"""``varus logan``: pre-screen SRA runs with Logan contig assemblies.

Logan (Chikhi et al. 2024) provides per-run assemblies of nearly every
public SRA run. It is rebuilt in full at discrete time points rather than
incrementally, so whether a run is covered is decided per run by an HTTP
``HEAD`` (:func:`check_availability`). Instead of blindly downloading reads
from hundreds of runs, ``varus logan`` downloads each candidate run's contig
FASTA (2-5 MB zstd each), spliced-aligns the contigs to the genome with
minimap2, and scores every run by the *breadth* of genome tiles it touches.
Foreign / contaminated / empty runs light up only a small fraction of the
tiles a genuine run does and are rejected; the remaining runs are ranked by
a greedy complementary-coverage selection (:mod:`varus.select`).

Outputs (all under ``<outdir>/logan/`` unless noted)::

    <outdir>/Runlist.logan.tsv   ranked runlist for the ``varus run`` stage
    LoganRanking.tsv             one row per candidate (all statuses)
    logan_introns.gff            stranded introns of the accepted runs
    logan.splice_sites           HISAT2 --known-splicesite-infile
    logan.junc.bed               minimap2 --junc-bed
    logan_tiles.tsv.gz           per-run tile weights of accepted runs
    logan_summary.json           parameters, counts, timings, coverage curve
    availability.tsv             cached HEAD results (resume)
    runs/<acc>.json + .tiles.npz per-run statistics (resume)

Conventions shared with the rest of VARUS
-----------------------------------------
* Tile index = ``(reference_start + 1) // tile_size`` (see :mod:`varus.tiles`).
* Introns are 1-based inclusive ``(chrom, start, end, strand)``; strand is
  assigned afterwards with :class:`varus.strand.StrandAssigner`.
* Contig read names are rewritten to ``<ACC>_<i>`` where ``i`` is the 0-based
  position of the contig in the file. The per-contig abundance array
  ``ka[i]`` uses the same index, which is why the original Logan counter is
  *not* preserved: Logan numbers contigs ``ACC_0, ACC_1, ...`` in file order
  anyway, so for well-formed files the two coincide.
"""

from __future__ import annotations

import concurrent.futures as cf
import gzip
import hashlib
import json
import logging
import math
import multiprocessing as mp
import os
import random
import shutil
import subprocess
import time
import urllib.error
import urllib.request
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from varus.align import (align_contigs_minimap2, minimap2_index_parts, remove_split_files,
                         reserve_threads, start_contig_alignment)
from varus.index import build_minimap2_index, check_sequence_lengths
from varus.introns import IntronCounts, IntronKey, iter_introns, write_introns_gff
from varus.runlist import RunRecord, write_runlist
from varus.strand import StrandAssigner, write_hisat2_splice_sites, write_minimap2_junc_bed
from varus.tiles import Tile

log = logging.getLogger(__name__)

LOGAN_CONTIGS_URL = "https://s3.amazonaws.com/logan-pub/c/{acc}/{acc}.contigs.fa.zst"
LOGAN_UNITIGS_URL = "https://s3.amazonaws.com/logan-pub/u/{acc}/{acc}.unitigs.fa.zst"


# opener(request, timeout) -> response; response must support .read(n) (for
# downloads), .getcode()/.status (for HEAD) and .close(). Injected in tests.
Opener = Callable[[urllib.request.Request, float], object]


def _default_opener(request: urllib.request.Request, timeout: float):
    return urllib.request.urlopen(request, timeout=timeout)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


@dataclass
class LoganConfig:
    """All knobs of the ``varus logan`` stage (mirrors the CLI options)."""

    genome: Path
    runlist: Path
    outdir: Path
    mmi: Optional[Path] = None
    threads: int = 4
    # Logan's S3 gave ~8 MB/s in total whether 8 or 16 connections were
    # open (Drosophila, 2026-09-24), so 8 is enough.
    download_workers: int = 8
    max_candidates: int = 500
    # Runs per minimap2 call; each call reloads the genome index (~2 s).
    chunk_runs: int = 25
    # minimap2's SAM is piped into this many scanner processes, so a chunk
    # is scanned while it is being aligned (0 = scan in the main process).
    # Only with --logan-bam are chunk BAMs written and scanned from disk.
    scan_workers: int = 2
    # minimap2 processes per chunk, each on a share of the chunk's runs and
    # threads and each piping into its own scanner. One scanner process parses
    # SAM slower than a 45-thread minimap2 writes it (25 Drosophila runs:
    # 38.7 s streamed vs 28.4 s alignment alone; 3 groups: 30.8 s). Results
    # are identical (minimap2 aligns every contig independently; a run is
    # never split). Every group loads its own copy of the index, so the
    # count is lowered to what fits in memory (max_groups_for_memory).
    # None = auto_align_groups(threads) (48 threads -> 3, <= 26 -> 1).
    align_groups: Optional[int] = None
    max_intron: int = 20_000
    min_contigs: int = 100
    min_tiles_frac: float = 0.10
    max_divergence: float = 0.05
    ka_cap: float = 50.0
    tile_size: int = 5000
    select_top: int = 50
    batch_size: int = 50_000
    prior_batches: float = 1.0
    keep_contigs: bool = False
    write_bam: bool = False
    seed: int = 1
    timeout: float = 60.0
    retries: int = 4

    @property
    def logan_dir(self) -> Path:
        return self.outdir / "logan"


# ---------------------------------------------------------------------------
# Header parsing and candidate sampling
# ---------------------------------------------------------------------------


def parse_logan_header(line: str) -> Tuple[str, float]:
    """Parse a Logan FASTA header line.

    ``>SRR10620183_0 ka:f:431.586   L:-:10472:-  `` -> ``("SRR10620183_0", 431.586)``.
    Both ``ka:f:`` and ``km:f:`` are accepted as the abundance tag. A missing
    or unparseable abundance yields ``math.nan``. Trailing whitespace and
    optional ``L:`` link fields are ignored.
    """
    s = line.strip()
    if s.startswith(">"):
        s = s[1:]
    parts = s.split()
    if not parts:
        return "", math.nan
    name = parts[0]
    ka = math.nan
    for tok in parts[1:]:
        if tok.startswith("ka:f:") or tok.startswith("km:f:"):
            try:
                ka = float(tok[5:])
            except ValueError:
                ka = math.nan
            break
    return name, ka


def auto_align_groups(threads: int) -> int:
    """Default ``--align-groups``: one minimap2 process per ~15 of its threads.

    One scanner keeps up with a 15-thread minimap2 but not with a 45-thread
    one (48 threads: 3 groups, the benchmarked value). At most 4, because
    beyond that the contig download (~8 MB/s on brain) is the limit and every
    group loads its own copy of the index.
    """
    return max(1, min(4, int((int(threads) - 4) / 15 + 0.5)))


def _read_int(path: Path) -> Optional[int]:
    try:
        v = path.read_text().strip()
    except OSError:
        return None
    return int(v) if v.isdigit() else None


def _cgroup_headroom(d: Path, limit_name: str, usage_name: str,
                     file_keys: Tuple[str, str]) -> Optional[int]:
    """Memory still free under the limit set in cgroup directory ``d``: the
    limit minus the usage, where the page cache (``file_keys`` in
    ``memory.stat``) counts as free because the kernel reclaims it first.
    None if ``d`` sets no limit."""
    limit = _read_int(d / limit_name)
    if limit is None or limit >= (1 << 60):
        return None
    usage = _read_int(d / usage_name)
    if usage is None:
        return limit
    cache = 0
    try:
        for line in (d / "memory.stat").read_text().splitlines():
            k, _, v = line.partition(" ")
            if k in file_keys and v.strip().isdigit():
                cache += int(v)
    except OSError:
        pass
    return max(0, min(limit, limit - usage + cache))


def _process_tree_rss() -> int:
    """Resident memory of this process and all its descendants, in bytes."""
    page = os.sysconf("SC_PAGE_SIZE")
    children: Dict[int, List[int]] = {}
    rss: Dict[int, int] = {}
    for d in Path("/proc").iterdir():
        if not d.name.isdigit():
            continue
        try:
            stat = (d / "stat").read_text()
        except OSError:
            continue
        fields = stat[stat.rindex(")") + 2:].split()   # after "pid (comm) "
        pid, ppid = int(d.name), int(fields[1])
        children.setdefault(ppid, []).append(pid)
        rss[pid] = int(fields[21]) * page
    total, todo = 0, [os.getpid()]
    while todo:
        pid = todo.pop()
        total += rss.get(pid, 0)
        todo.extend(children.get(pid, ()))
    return total


def _slurm_headroom() -> Optional[int]:
    """The SLURM job's memory allocation (--mem or --mem-per-cpu, MiB) minus
    what this process tree holds. The job's cgroup is not visible inside a
    Singularity container (/sys/fs/cgroup is empty there; brain, 2026-09-26),
    but SLURM's environment is."""
    env = os.environ
    try:
        if env.get("SLURM_MEM_PER_NODE"):
            limit = int(env["SLURM_MEM_PER_NODE"]) * 2**20
        elif env.get("SLURM_MEM_PER_CPU"):
            cpus = int(env.get("SLURM_CPUS_ON_NODE") or env.get("SLURM_CPUS_PER_TASK") or 1)
            limit = int(env["SLURM_MEM_PER_CPU"]) * 2**20 * cpus
        else:
            return None
    except ValueError:
        return None
    try:
        used = _process_tree_rss()
    except (OSError, ValueError, IndexError):
        used = 0
    return max(0, limit - used)


def available_memory_bytes() -> Optional[int]:
    """Memory this process may still use: the smallest of the node's
    MemAvailable, the headroom under the tightest cgroup limit and the SLURM
    allocation minus this process tree (the cgroup is invisible inside a
    container). Measured now, so it drops while other processes of the job
    or the node hold memory."""
    limits: List[int] = []
    h = _slurm_headroom()
    if h is not None:
        limits.append(h)
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                limits.append(int(line.split()[1]) * 1024)
    except (OSError, ValueError, IndexError):
        pass
    try:
        cgroups = Path("/proc/self/cgroup").read_text().splitlines()
    except OSError:
        cgroups = []
    for line in cgroups:
        parts = line.split(":", 2)
        if len(parts) != 3:
            continue
        _, ctrls, path = parts
        if ctrls == "":                                   # cgroup v2
            base = Path("/sys/fs/cgroup")
            names = ("memory.max", "memory.current", ("active_file", "inactive_file"))
        elif "memory" in ctrls.split(","):                # cgroup v1
            base = Path("/sys/fs/cgroup/memory")
            names = ("memory.limit_in_bytes", "memory.usage_in_bytes",
                     ("total_active_file", "total_inactive_file"))
        else:
            continue
        rel = Path(path.strip("/"))
        for d in [rel, *rel.parents]:                     # limits may sit on a parent
            h = _cgroup_headroom(base / d, *names)
            if h is not None:
                limits.append(h)
    return min(limits) if limits else None


# Memory of one minimap2 process besides its index: the query batch and
# per-thread buffers.
_GROUP_OVERHEAD = 2 * 2**30


def index_memory_bytes(mmi: Path) -> Optional[int]:
    """Index memory of one minimap2 process: the largest part of ``mmi``
    (minimap2 loads one part at a time), or the file size if the parts
    cannot be read. None if ``mmi`` does not exist."""
    parts = minimap2_index_parts(Path(mmi))
    if parts:
        return max(parts)
    try:
        return Path(mmi).stat().st_size
    except OSError:
        return None


def max_groups_for_memory(mmi: Path, avail: Optional[int] = None) -> int:
    """How many minimap2 processes fit in memory, each with its own copy of
    ``mmi`` (its largest part) plus ~2 GB of query batch and buffers, within
    60 % of what is available. Unknown memory: one process for indexes over 4 GB."""
    idx = index_memory_bytes(mmi)
    if idx is None:
        return 1 << 10
    if avail is None:
        avail = available_memory_bytes()
    if avail is None:
        return 1 if idx > 4 * 2**30 else 1 << 10
    return max(1, int(0.6 * avail // (idx + _GROUP_OVERHEAD)))


# minimap2 -x splice index per genome base and the peak of building it in
# one part, measured on wheat (14.6 Gbp: 51.4 GB index, 88.5 GB build peak,
# 2026-09-25) and mouse (4.0 bytes/base). minimap2 splits every 8 Gbp.
_INDEX_BYTES_PER_BASE = 4.0
_INDEX_BUILD_PEAK = 2.0            # x index size; wheat 1.72
_MM2_PART_BASES = 8_000_000_000


def single_part_bases(genome: Path, avail: Optional[int] = None) -> Optional[int]:
    """minimap2 ``-I`` that keeps the index of ``genome`` in one part, or None
    for minimap2's default.

    A genome over 8 Gbp is split into parts by default, and a split index
    aligns every contig once per part: wheat took 5.9 h split against 2.6 h in
    one part, with the same 50 runs selected (2026-09-26). One part is used
    when building it (~8 bytes per base) fits in 90 % of the available
    memory. The base count is bounded by the file size (x5 if gzipped), which
    is also the ``-I`` passed.
    """
    try:
        size = Path(genome).stat().st_size
    except OSError:
        return None
    bases = size * 5 if str(genome).endswith(".gz") else size
    if bases <= _MM2_PART_BASES:
        return None                              # one part anyway
    need = _INDEX_BUILD_PEAK * _INDEX_BYTES_PER_BASE * bases
    if avail is None:
        avail = available_memory_bytes()
    if avail is None or need > 0.9 * avail:
        log.info("minimap2 index of %s in parts of 8 Gbp: one part needs ~%.0f GB to build, "
                 "%s available", genome, need / 1e9,
                 "unknown" if avail is None else f"{avail / 1e9:.0f} GB")
        return None
    log.info("minimap2 index of %s in one part (-I %d): ~%.0f GB to build, %.0f GB available",
             genome, bases, need / 1e9, avail / 1e9)
    return bases


def warn_if_index_does_not_fit(mmi: Path, avail: Optional[int] = None) -> bool:
    """Log a warning if even one minimap2 process (index + overhead) exceeds
    the available memory. Returns True if it fits or memory is unknown."""
    idx = index_memory_bytes(mmi)
    if idx is None:
        return True
    need = idx + _GROUP_OVERHEAD
    if avail is None:
        avail = available_memory_bytes()
    if avail is None or need <= avail:
        return True
    log.warning("minimap2 needs ~%.1f GB for the index %s but only %.1f GB of memory is "
                "available; expect swapping or an out-of-memory kill. Request more memory "
                "(e.g. SLURM --mem).", need / 2**30, mmi, avail / 2**30)
    return False


def _sha1(s: str) -> str:
    return hashlib.sha1(s.encode("utf-8")).hexdigest()


def sample_candidates(records: List[RunRecord], max_candidates: int) -> List[RunRecord]:
    """Deterministic, project-balanced ordering of candidate runs.

    Bioprojects are ordered by ``sha1(bioproject)`` (a run without a
    bioproject forms its own bucket keyed by its accession); within a
    bioproject runs are ordered by ``sha1(accession)``; the buckets are then
    visited round-robin so the first ``max_candidates`` runs span as many
    projects as possible. ``max_candidates <= 0`` returns the full ordering.
    Pure function; the same input always yields the same output.
    """
    buckets: Dict[str, List[RunRecord]] = {}
    for r in records:
        key = r.bioproject if r.bioproject else f"\x00{r.accession}"
        buckets.setdefault(key, []).append(r)
    ordered_keys = sorted(buckets, key=_sha1)
    for k in ordered_keys:
        buckets[k].sort(key=lambda r: _sha1(r.accession))
    out: List[RunRecord] = []
    depth = 0
    while True:
        added = False
        for k in ordered_keys:
            b = buckets[k]
            if depth < len(b):
                out.append(b[depth])
                added = True
        if not added:
            break
        depth += 1
    if max_candidates and max_candidates > 0:
        out = out[:max_candidates]
    return out


def read_runlist(path: Path) -> List[RunRecord]:
    """Parse ``Runlist.tsv`` into :class:`RunRecord` objects.

    Same rules as :func:`varus.controller.load_runs`: lines starting with
    ``@`` are the header, 6-8 tab-separated columns, malformed lines are
    skipped with a warning, colorspace runs are dropped.
    """
    records: List[RunRecord] = []
    with Path(path).open(encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("@"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            try:
                spots = int(parts[1])
                bases = int(parts[2])
                avg_len = float(parts[3])
                paired = bool(int(parts[4]))
                colorspace = bool(int(parts[5]))
            except (ValueError, IndexError):
                log.warning("Skipping malformed runlist line: %s", line.rstrip())
                continue
            if colorspace:
                continue
            records.append(RunRecord(
                accession=parts[0].strip(),
                total_spots=spots,
                total_bases=bases,
                avg_len=avg_len,
                paired=paired,
                colorspace=colorspace,
                platform=parts[6].strip() if len(parts) >= 7 else "",
                bioproject=parts[7].strip() if len(parts) >= 8 else "",
            ))
    log.info("Read %d runs from %s", len(records), path)
    return records


# ---------------------------------------------------------------------------
# Availability (HTTP HEAD)
# ---------------------------------------------------------------------------


def _backoff(attempt: int, base: float = 1.0, cap: float = 30.0) -> float:
    return min(cap, base * (2 ** attempt)) * (0.5 + random.random())


def head_available(
    acc: str,
    *,
    timeout: float,
    retries: int,
    url_template: str = LOGAN_CONTIGS_URL,
    opener: Optional[Opener] = None,
) -> Tuple[str, int]:
    """HEAD the Logan object for ``acc``; return ``(status, http_code)``.

    ``status`` is ``"available"`` (2xx), ``"absent"`` (404 or 403 - S3 answers
    403 for missing keys when listing is denied) or ``"error"`` (5xx, network
    errors or timeouts after ``retries`` attempts with exponential backoff).
    """
    opener = opener or _default_opener
    url = url_template.format(acc=acc)
    last_code = 0
    attempts = max(1, int(retries) + 1)
    for attempt in range(attempts):
        req = urllib.request.Request(url, method="HEAD")
        try:
            resp = opener(req, timeout)
            try:
                code = getattr(resp, "status", None)
                if code is None:
                    code = resp.getcode()
            finally:
                close = getattr(resp, "close", None)
                if close is not None:
                    close()
            if 200 <= int(code) < 300:
                return "available", int(code)
            last_code = int(code)
        except urllib.error.HTTPError as e:
            last_code = int(e.code)
            if e.code in (404, 403):
                return "absent", int(e.code)
            if e.code < 500:
                log.warning("HEAD %s: HTTP %d", url, e.code)
                return "error", int(e.code)
        except (urllib.error.URLError, OSError, TimeoutError) as e:
            log.debug("HEAD %s attempt %d failed: %s", url, attempt + 1, e)
            last_code = 0
        if attempt < attempts - 1:
            time.sleep(_backoff(attempt))
    return "error", last_code


def _read_availability_cache(path: Path) -> Dict[str, Tuple[str, int]]:
    out: Dict[str, Tuple[str, int]] = {}
    if not path.is_file():
        return out
    with path.open(encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] == "acc":
                continue
            try:
                out[parts[0]] = (parts[1], int(parts[2] or 0))
            except ValueError:
                continue
    return out


def check_availability(
    accs: List[str],
    cache_path: Path,
    *,
    workers: int,
    timeout: float,
    retries: int,
    url_template: str = LOGAN_CONTIGS_URL,
    opener: Optional[Opener] = None,
) -> Dict[str, Tuple[str, int]]:
    """Check which accessions have a Logan contig file, with an on-disk cache.

    The cache (``acc\\tstatus\\thttp\\tchecked_at``) is read first; only
    accessions missing from it or previously recorded as ``error`` are
    queried (in parallel), and the results are appended. Returns the mapping
    for all requested accessions.
    """
    cache = _read_availability_cache(cache_path)
    todo = [a for a in accs if a not in cache or cache[a][0] == "error"]
    if todo:
        log.info("HEAD-checking %d/%d accessions on Logan S3 (%d cached)",
                 len(todo), len(accs), len(accs) - len(todo))
        results: Dict[str, Tuple[str, int]] = {}
        with cf.ThreadPoolExecutor(max_workers=max(1, workers)) as ex:
            futs = {
                ex.submit(head_available, a, timeout=timeout, retries=retries,
                          url_template=url_template, opener=opener): a
                for a in todo
            }
            for fut in cf.as_completed(futs):
                a = futs[fut]
                try:
                    results[a] = fut.result()
                except Exception as e:  # pragma: no cover - defensive
                    log.warning("HEAD %s raised %s", a, e)
                    results[a] = ("error", 0)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        new_file = not cache_path.exists()
        now = datetime.now(timezone.utc).isoformat(timespec="seconds")
        with cache_path.open("a", encoding="utf-8") as f:
            if new_file:
                f.write("acc\tstatus\thttp\tchecked_at\n")
            for a in todo:  # keep request order for readability
                st, code = results[a]
                f.write(f"{a}\t{st}\t{code}\t{now}\n")
        cache.update(results)
    return {a: cache[a] for a in accs}


# ---------------------------------------------------------------------------
# Contig download
# ---------------------------------------------------------------------------


class LoganDownloadError(RuntimeError):
    """Download failed after all retries (network / server / decode error)."""

    def __init__(self, msg: str, http: int = 0) -> None:
        super().__init__(msg)
        self.http = http


class LoganAbsentError(LoganDownloadError):
    """The object does not exist on Logan S3 (404/403)."""


@dataclass
class LoganContigs:
    acc: str
    fasta: Path
    ka: np.ndarray  # float32, one entry per contig (NaN if unparseable)
    n_contigs: int
    total_bp: int
    bytes_downloaded: int
    seconds: float


class _CountingReader:
    """File-like wrapper that counts the bytes read from ``raw``."""

    def __init__(self, raw) -> None:
        self.raw = raw
        self.n = 0

    def read(self, size: int = -1) -> bytes:
        b = self.raw.read(size)
        self.n += len(b)
        return b


def _iter_decompressed(counting: _CountingReader, acc: str):
    """Yield decompressed chunks of a zstd HTTP body.

    Uses ``zstandard`` when importable (handles concatenated frames and
    raises :class:`LoganDownloadError` when the stream ends inside a frame,
    i.e. the download was truncated - the plain ``stream_reader`` would
    silently return what it has). Falls back to a ``zstd -dc`` subprocess
    when the package is missing.
    """
    try:
        import zstandard  # optional extra "logan"
    except ImportError:  # pragma: no cover - exercised only without zstandard
        zstandard = None
    if zstandard is not None:
        dctx = zstandard.ZstdDecompressor()
        dobj = dctx.decompressobj()
        fed = 0
        while True:
            chunk = counting.read(1 << 20)
            if not chunk:
                break
            while chunk:
                out = dobj.decompress(chunk)
                fed += len(chunk)
                if out:
                    yield out
                if dobj.eof:
                    chunk = dobj.unused_data
                    dobj = dctx.decompressobj()
                    fed = 0
                else:
                    chunk = b""
        if fed > 0 and not dobj.eof:
            raise LoganDownloadError(f"{acc}: truncated zstd stream")
        return

    if shutil.which("zstd") is None:
        raise LoganDownloadError("neither the 'zstandard' package nor a 'zstd' binary is available")
    import threading

    proc = subprocess.Popen(["zstd", "-dc"], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    feed_err: List[BaseException] = []

    def _feed() -> None:
        try:
            while True:
                chunk = counting.read(1 << 20)
                if not chunk:
                    break
                proc.stdin.write(chunk)
        except BaseException as e:  # propagated after the consumer finishes
            feed_err.append(e)
        finally:
            try:
                proc.stdin.close()
            except OSError:
                pass

    t = threading.Thread(target=_feed, daemon=True)
    t.start()
    try:
        while True:
            out = proc.stdout.read(1 << 20)
            if not out:
                break
            yield out
        t.join()
        rc = proc.wait()
        if feed_err:
            raise feed_err[0]
        if rc != 0:
            raise LoganDownloadError(f"{acc}: zstd -dc exited with status {rc}")
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()


def _load_contigs_sidecar(acc: str, fasta: Path, ka_path: Path) -> LoganContigs:
    ka = np.load(ka_path).astype(np.float32, copy=False)
    total_bp = 0
    with fasta.open("rb") as f:
        for line in f:
            if not line.startswith(b">"):
                total_bp += len(line.rstrip(b"\r\n"))
    return LoganContigs(acc=acc, fasta=fasta, ka=ka, n_contigs=int(ka.size),
                        total_bp=total_bp, bytes_downloaded=0, seconds=0.0)


def download_contigs(
    acc: str,
    dest_dir: Path,
    *,
    timeout: float,
    retries: int,
    url_template: str = LOGAN_CONTIGS_URL,
    opener: Optional[Opener] = None,
) -> LoganContigs:
    """Stream ``<acc>.contigs.fa.zst`` from Logan and write a plain FASTA.

    The FASTA is written to ``<dest_dir>/<acc>.contigs.fa`` with every header
    rewritten to ``>ACC_<i>`` (``i`` = 0-based running index; see the module
    docstring for why the original Logan counter is not kept). The per-contig
    abundance (``ka:f:`` / ``km:f:``) is saved next to it as
    ``<acc>.ka.npy``; when both files already exist the download is skipped
    (resume). Output is written to a ``.part`` file and renamed on success.

    Raises :class:`LoganAbsentError` on 404/403 and
    :class:`LoganDownloadError` after ``retries`` failed attempts.
    """
    opener = opener or _default_opener
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    fasta = dest_dir / f"{acc}.contigs.fa"
    ka_path = dest_dir / f"{acc}.ka.npy"
    if fasta.is_file() and ka_path.is_file():
        log.debug("%s: reusing existing contigs %s", acc, fasta)
        return _load_contigs_sidecar(acc, fasta, ka_path)

    url = url_template.format(acc=acc)
    part = dest_dir / f"{acc}.contigs.fa.part"
    attempts = max(1, int(retries) + 1)
    last_err: Optional[Exception] = None
    for attempt in range(attempts):
        t0 = time.monotonic()
        resp = None
        try:
            req = urllib.request.Request(url, method="GET")
            resp = opener(req, timeout)
            counting = _CountingReader(resp)
            kas: List[float] = []
            n_contigs = 0
            total_bp = 0
            warned = False
            buf = b""

            def _emit(line: bytes) -> None:
                nonlocal n_contigs, total_bp, warned
                if line.startswith(b">"):
                    _, ka = parse_logan_header(line.decode("utf-8", "replace"))
                    if math.isnan(ka) and not warned:
                        log.warning("%s: contig header without parseable abundance: %s",
                                    acc, line[:60].decode("utf-8", "replace"))
                        warned = True
                    kas.append(ka)
                    out.write(b">%s_%d\n" % (acc.encode(), n_contigs))
                    n_contigs += 1
                else:
                    seq = line.rstrip(b"\r")
                    total_bp += len(seq)
                    out.write(seq + b"\n")

            with part.open("wb") as out:
                for chunk in _iter_decompressed(counting, acc):
                    buf += chunk
                    lines = buf.split(b"\n")
                    buf = lines.pop()
                    for line in lines:
                        _emit(line)
                if buf:
                    _emit(buf)
            ka_arr = np.asarray(kas, dtype=np.float32)
            np.save(ka_path, ka_arr)
            os.replace(part, fasta)
            secs = time.monotonic() - t0
            mb = counting.n / 1e6
            log.info("%s: %d contigs, %.1f Mbp, %.1f MB in %.1f s (%.1f MB/s)",
                     acc, n_contigs, total_bp / 1e6, mb, secs, mb / secs if secs > 0 else 0.0)
            return LoganContigs(acc=acc, fasta=fasta, ka=ka_arr, n_contigs=n_contigs,
                                total_bp=total_bp, bytes_downloaded=counting.n, seconds=secs)
        except urllib.error.HTTPError as e:
            if e.code in (404, 403):
                raise LoganAbsentError(f"{url}: HTTP {e.code}", http=e.code) from e
            last_err = e
            log.warning("%s: HTTP %d on attempt %d/%d", acc, e.code, attempt + 1, attempts)
        except Exception as e:  # URLError, OSError, timeouts, zstd errors
            last_err = e
            log.warning("%s: download attempt %d/%d failed: %s", acc, attempt + 1, attempts, e)
        finally:
            if resp is not None:
                close = getattr(resp, "close", None)
                if close is not None:
                    try:
                        close()
                    except Exception:  # pragma: no cover
                        pass
            if part.exists():
                try:
                    part.unlink()
                except OSError:  # pragma: no cover
                    pass
        if attempt < attempts - 1:
            time.sleep(_backoff(attempt))
    http = getattr(last_err, "code", 0) or 0
    raise LoganDownloadError(f"{acc}: download failed after {attempts} attempts: {last_err}", http=int(http))


# ---------------------------------------------------------------------------
# Per-run statistics
# ---------------------------------------------------------------------------


@dataclass
class LoganRunStats:
    acc: str
    n_contigs: int = 0
    total_bp: int = 0
    n_aligned: int = 0
    aligned_bp: int = 0
    mapped_pct: float = 0.0
    # Abundance x length mass of contigs (all / mapped): approximates the
    # fraction of the run's *reads* that come from the target genome, which
    # the plain contig-count fraction does not (many unmapped contigs are
    # low-abundance junk). yield_pct = 100 * mass_mapped / mass_total.
    mass_total: float = 0.0
    mass_mapped: float = 0.0
    yield_pct: float = 0.0
    n_spliced: int = 0
    n_tiles: int = 0
    tile_mass: float = 0.0
    # Aligned-bp weighted median of minimap2's gap-compressed divergence
    # (``de`` tag) over primary contig alignments; NaN = unknown (nothing
    # aligned, or stats cached before this field existed). Contigs from the
    # target species sit near 0; other species of the same genus at 0.1-0.2,
    # where their contigs still cover the genome but HISAT2 maps < 5 % of reads.
    divergence: float = math.nan
    tiles: Dict[Tile, float] = field(default_factory=dict)
    introns: Dict[IntronKey, float] = field(default_factory=dict)  # strand "."
    status: str = "pending"  # accepted|rejected|too_few_contigs|absent|error|unsampled
    rank: int = 0
    gain: float = 0.0
    http: int = 0
    bioproject: str = ""


def _contig_weight(ka: float, cap: float) -> float:
    if ka is None or math.isnan(ka):
        return 1.0
    return float(min(ka, cap))


def scan_chunk_bam(
    bam,
    ka_by_acc: Dict[str, np.ndarray],
    n_contigs_by_acc: Dict[str, int],
    total_bp_by_acc: Dict[str, int],
    *,
    tile_size: int,
    ka_cap: float,
) -> Dict[str, LoganRunStats]:
    """One pysam pass over a chunk's alignments; per-run tile and intron weights.

    ``bam`` is a path to a BAM/SAM file or a binary file object (e.g. the
    read end of a pipe fed by minimap2; the format is auto-detected and the
    stream is read once, in order, so no index or sorting is needed).

    Only primary alignments are counted (unmapped, secondary and
    supplementary records are skipped). The run is recovered from the query
    name ``<ACC>_<i>``. Each aligned contig adds
    ``min(ka, ka_cap) * max(1, aligned_ref_len / 150)`` to its tile (NaN
    abundance counts as 1), i.e. roughly the number of 150-bp reads it
    stands for; introns (CIGAR ``N``) are weighted with ``min(ka, ka_cap)``.
    Every accession in ``ka_by_acc`` gets a stats entry even if none of its
    contigs aligned.
    """
    if tile_size <= 0:
        raise ValueError("tile_size must be > 0")
    import pysam  # optional extra "align"

    stats: Dict[str, LoganRunStats] = {}
    for acc in ka_by_acc:
        stats[acc] = LoganRunStats(
            acc=acc,
            n_contigs=int(n_contigs_by_acc.get(acc, len(ka_by_acc[acc]))),
            total_bp=int(total_bp_by_acc.get(acc, 0)),
        )
    unknown: set = set()
    div: Dict[str, List[Tuple[float, int]]] = {acc: [] for acc in stats}
    is_path = isinstance(bam, (str, Path))
    with pysam.AlignmentFile(str(bam) if is_path else bam, "rb" if is_path else "r") as fh:
        for read in fh.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            qname = read.query_name or ""
            acc, _, idx_s = qname.rpartition("_")
            if acc not in stats:
                if acc not in unknown:
                    log.warning("BAM %s: read %s from unknown run; ignored", bam, qname)
                    unknown.add(acc)
                continue
            st = stats[acc]
            ka_arr = ka_by_acc[acc]
            try:
                idx = int(idx_s)
                ka = float(ka_arr[idx]) if 0 <= idx < len(ka_arr) else math.nan
            except ValueError:
                ka = math.nan
            # Read mass proxy: raw abundance (uncapped, NaN -> 1) x contig length,
            # accumulated for every contig so mapped / total is a yield estimate.
            ka_raw = 1.0 if (ka != ka) else max(ka, 0.0)
            qlen = read.query_length or read.infer_read_length() or 0
            mass = ka_raw * float(qlen)
            st.mass_total += mass
            if read.is_unmapped:
                continue
            st.mass_mapped += mass
            wk = _contig_weight(ka, ka_cap)
            ref_len = read.reference_end - read.reference_start if read.reference_end is not None else 0
            w = wk * max(1.0, ref_len / 150.0)
            tile: Tile = (read.reference_name, (read.reference_start + 1) // tile_size)
            st.tiles[tile] = st.tiles.get(tile, 0.0) + w
            st.n_aligned += 1
            alen = int(read.query_alignment_length)
            st.aligned_bp += alen
            d = _read_divergence(read)
            if d is not None and alen > 0:
                div[acc].append((d, alen))
            spliced = False
            for chrom, s, e in iter_introns(read):
                key: IntronKey = (chrom, s, e, ".")
                st.introns[key] = st.introns.get(key, 0.0) + wk
                spliced = True
            if spliced:
                st.n_spliced += 1
    for st in stats.values():
        st.n_tiles = len(st.tiles)
        st.tile_mass = float(sum(st.tiles.values()))
        st.mapped_pct = 100.0 * st.n_aligned / st.n_contigs if st.n_contigs > 0 else 0.0
        st.yield_pct = 100.0 * st.mass_mapped / st.mass_total if st.mass_total > 0 else 0.0
        st.divergence = _weighted_median(div[st.acc])
    return stats


def _scan_chunk_timed(*args, **kwargs) -> Tuple[Dict[str, LoganRunStats], float]:
    """:func:`scan_chunk_bam` plus its own wall time; runs in a scan worker."""
    t0 = time.monotonic()
    res = scan_chunk_bam(*args, **kwargs)
    return res, time.monotonic() - t0


def _scan_stream_timed(conn, *args, **kwargs) -> Tuple[Dict[str, LoganRunStats], float]:
    """Scan the SAM that minimap2 writes into pipe ``conn`` (read end).

    ``conn`` is a :class:`multiprocessing.connection.Connection` wrapping
    the raw pipe descriptor; it is used only to carry the descriptor into a
    spawned scan worker (the parent keeps its copy open until the result is
    collected, because the pool pickles the task lazily). The wall time
    returned includes waiting for minimap2's output.
    """
    t0 = time.monotonic()
    fh = os.fdopen(os.dup(conn.fileno()), "rb")
    conn.close()
    try:
        res = scan_chunk_bam(fh, *args, **kwargs)
    finally:
        fh.close()
    return res, time.monotonic() - t0


def _minimap2_summary(err: Path) -> str:
    """minimap2's last stderr line (``Real time ...; CPU ...``), or ''."""
    try:
        lines = err.read_text(errors="replace").splitlines()
    except OSError:
        return ""
    return next((l.split("] ", 1)[-1] for l in reversed(lines) if "Real time" in l), "")


def _read_divergence(read) -> Optional[float]:
    """minimap2 ``de`` tag; falls back to NM / aligned length."""
    if read.has_tag("de"):
        return float(read.get_tag("de"))
    if read.has_tag("NM") and read.query_alignment_length:
        return float(read.get_tag("NM")) / float(read.query_alignment_length)
    return None


def _weighted_median(pairs: List[Tuple[float, int]]) -> float:
    """Median of values weighted by ``w``; NaN for an empty list."""
    if not pairs:
        return math.nan
    pairs = sorted(pairs)
    half = sum(w for _, w in pairs) / 2.0
    acc = 0
    for v, w in pairs:
        acc += w
        if acc >= half:
            return float(v)
    return float(pairs[-1][0])


def apply_gate(
    stats: Dict[str, LoganRunStats],
    *,
    min_contigs: int,
    min_tiles_frac: float,
    max_divergence: float = 0.05,
) -> None:
    """Set ``status`` of every scanned run in place.

    * ``too_few_contigs`` if ``n_contigs < min_contigs``,
    * ``rejected`` if the median contig divergence exceeds ``max_divergence``
      (another species: its contigs cover the genome under minimap2, but its
      reads fail HISAT2; unknown divergence passes),
    * otherwise ``accepted`` if ``n_tiles >= min_tiles_frac * max_tiles``
      (``max_tiles`` = best ``n_tiles`` among runs with enough contigs that
      pass the divergence check) and ``n_tiles > 0``, else ``rejected``.

    Runs whose status is ``absent`` / ``error`` / ``unsampled`` are untouched.
    """
    scannable = [s for s in stats.values()
                 if s.status in ("pending", "accepted", "rejected", "too_few_contigs")]
    def _divergent(s: LoganRunStats) -> bool:
        return max_divergence > 0 and s.divergence == s.divergence and s.divergence > max_divergence

    eligible = [s for s in scannable if s.n_contigs >= min_contigs and not _divergent(s)]
    max_tiles = max((s.n_tiles for s in eligible), default=0)
    thr = min_tiles_frac * max_tiles
    n_acc = n_rej = n_few = n_div = 0
    for s in scannable:
        if s.n_contigs < min_contigs:
            s.status = "too_few_contigs"
            n_few += 1
        elif _divergent(s):
            s.status = "rejected"
            n_rej += 1
            n_div += 1
        elif s.n_tiles > 0 and s.n_tiles >= thr:
            s.status = "accepted"
            n_acc += 1
        else:
            s.status = "rejected"
            n_rej += 1
    log.info("Gate: max_tiles=%d threshold=%.1f max_divergence=%.3f -> %d accepted, "
             "%d rejected (%d divergent), %d too few contigs",
             max_tiles, thr, max_divergence, n_acc, n_rej, n_div, n_few)


# ---------------------------------------------------------------------------
# Persistence of per-run stats (resume)
# ---------------------------------------------------------------------------

_SCALAR_FIELDS = ("acc", "n_contigs", "total_bp", "n_aligned", "aligned_bp", "mapped_pct",
                  "mass_total", "mass_mapped", "yield_pct",
                  "n_spliced", "n_tiles", "tile_mass", "divergence", "status", "rank", "gain", "http",
                  "bioproject")


def save_run_stats(stats: LoganRunStats, runs_dir: Path) -> Tuple[Path, Path]:
    """Write ``<acc>.json`` (scalars) and ``<acc>.tiles.npz`` (arrays)."""
    runs_dir = Path(runs_dir)
    runs_dir.mkdir(parents=True, exist_ok=True)
    jpath = runs_dir / f"{stats.acc}.json"
    npath = runs_dir / f"{stats.acc}.tiles.npz"
    tiles = list(stats.tiles.items())
    introns = list(stats.introns.items())
    np.savez(
        npath,
        tile_chrom=np.asarray([t[0][0] for t in tiles], dtype=str),
        tile_idx=np.asarray([t[0][1] for t in tiles], dtype=np.int64),
        tile_w=np.asarray([t[1] for t in tiles], dtype=np.float32),
        intron_chrom=np.asarray([k[0] for k, _ in introns], dtype=str),
        intron_start=np.asarray([k[1] for k, _ in introns], dtype=np.int64),
        intron_end=np.asarray([k[2] for k, _ in introns], dtype=np.int64),
        intron_w=np.asarray([w for _, w in introns], dtype=np.float32),
    )
    # JSON has no NaN; None round-trips to NaN in load_run_stats
    scalars = {k: getattr(stats, k) for k in _SCALAR_FIELDS}
    if scalars["divergence"] != scalars["divergence"]:
        scalars["divergence"] = None
    tmp = jpath.with_suffix(".json.part")
    tmp.write_text(json.dumps(scalars, indent=1), encoding="utf-8")
    os.replace(tmp, jpath)
    return jpath, npath


def load_run_stats(acc: str, runs_dir: Path) -> Optional[LoganRunStats]:
    """Inverse of :func:`save_run_stats`; ``None`` if the JSON is missing."""
    runs_dir = Path(runs_dir)
    jpath = runs_dir / f"{acc}.json"
    npath = runs_dir / f"{acc}.tiles.npz"
    if not jpath.is_file():
        return None
    try:
        scalars = json.loads(jpath.read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        log.warning("Could not read %s (%s); will recompute", jpath, e)
        return None
    st = LoganRunStats(acc=acc)
    for k in _SCALAR_FIELDS:
        if k in scalars and scalars[k] is not None:
            setattr(st, k, scalars[k])
    if npath.is_file():
        with np.load(npath) as z:
            for chrom, idx, w in zip(z["tile_chrom"], z["tile_idx"], z["tile_w"]):
                st.tiles[(str(chrom), int(idx))] = float(w)
            for chrom, s, e, w in zip(z["intron_chrom"], z["intron_start"], z["intron_end"], z["intron_w"]):
                st.introns[(str(chrom), int(s), int(e), ".")] = float(w)
    return st


# ---------------------------------------------------------------------------
# Prior for the ``varus run`` stage
# ---------------------------------------------------------------------------


@dataclass
class LoganPrior:
    status: Dict[str, str]
    rank: Dict[str, int]
    tiles: Dict[str, Dict[Tile, float]]
    introns: Optional[Path]
    splice_sites: Optional[Path]
    junc_bed: Optional[Path]
    params: dict
    # yield_pct per accepted run (read-mass fraction on target, 0-100)
    yield_pct: Dict[str, float] = field(default_factory=dict)
    # status counts from logan_summary.json (accepted/rejected/...)
    counts: Dict[str, int] = field(default_factory=dict)

    @property
    def acceptance_rate(self) -> Optional[float]:
        """Fraction of *scanned* runs the gate accepted; None if unknown."""
        acc = self.counts.get("accepted", 0)
        scanned = acc + self.counts.get("rejected", 0) + self.counts.get("too_few_contigs", 0)
        return acc / scanned if scanned > 0 else None


def load_logan(logan_dir: Path) -> LoganPrior:
    """Read the outputs of :func:`run_logan`; missing files are tolerated."""
    logan_dir = Path(logan_dir)
    status: Dict[str, str] = {}
    rank: Dict[str, int] = {}
    yield_pct: Dict[str, float] = {}
    ranking = logan_dir / "LoganRanking.tsv"
    if ranking.is_file():
        with ranking.open(encoding="utf-8") as f:
            header = f.readline().rstrip("\n").split("\t")
            col = {name: i for i, name in enumerate(header)}

            def _get(parts, name):
                i = col.get(name)
                return parts[i] if i is not None and i < len(parts) else ""

            for line in f:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                acc = parts[col["acc"]]
                status[acc] = parts[col["status"]]
                r = _get(parts, "selected_rank")
                rank[acc] = int(r) if r else 0
                y = _get(parts, "yield_pct")
                if y:
                    try:
                        yield_pct[acc] = float(y)
                    except ValueError:
                        pass
    tiles: Dict[str, Dict[Tile, float]] = {}
    tpath = logan_dir / "logan_tiles.tsv.gz"
    if tpath.is_file():
        with gzip.open(tpath, "rt", encoding="utf-8") as f:
            for line in f:
                if not line.strip() or line.startswith("acc\t"):
                    continue
                acc, chrom, idx, w = line.rstrip("\n").split("\t")[:4]
                tiles.setdefault(acc, {})[(chrom, int(idx))] = float(w)
    params: dict = {}
    counts: Dict[str, int] = {}
    spath = logan_dir / "logan_summary.json"
    if spath.is_file():
        try:
            summary = json.loads(spath.read_text(encoding="utf-8"))
            params = summary.get("params", {}) or {}
            counts = {k: int(v) for k, v in (summary.get("counts", {}) or {}).items()}
        except ValueError:
            params = {}

    def _opt(name: str) -> Optional[Path]:
        p = logan_dir / name
        return p if p.is_file() else None

    return LoganPrior(
        status=status, rank=rank, tiles=tiles,
        introns=_opt("logan_introns.gff"),
        splice_sites=_opt("logan.splice_sites"),
        junc_bed=_opt("logan.junc.bed"),
        params=params,
        yield_pct=yield_pct,
        counts=counts,
    )


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def _fmt(x, nd: int = 2) -> str:
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    return str(x)


def _write_ranking(path: Path, order: List[RunRecord], stats: Dict[str, LoganRunStats],
                   curve: List[Tuple[str, float, float, int]]) -> None:
    """Write LoganRanking.tsv: selected runs first (by rank), then the other
    accepted runs by tile breadth, then rejected / too_few_contigs / absent /
    error / unsampled candidates."""
    cum_s = {acc: s for acc, _, s, _ in curve}
    cols = ["acc", "bioproject", "status", "http", "n_contigs", "contig_mb", "mapped_pct",
            "yield_pct", "divergence", "n_tiles", "tile_mass", "n_spliced", "n_introns", "selected_rank",
            "gain", "cumulative_S"]
    status_order = {s: i for i, s in enumerate(
        ("accepted", "rejected", "too_few_contigs", "absent", "error", "unsampled"))}

    def _key(rec: RunRecord):
        st = stats.get(rec.accession)
        if st is None:
            return (10**9, len(status_order), 0, rec.accession)
        return (st.rank or 10**9, status_order.get(st.status, len(status_order)),
                -st.n_tiles, rec.accession)

    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for rec in sorted(order, key=_key):
            st = stats.get(rec.accession)
            if st is None:
                f.write(f"{rec.accession}\t{rec.bioproject}\tunsampled\t" + "\t" * 13 + "\n")
                continue
            scanned = st.status in ("accepted", "rejected", "too_few_contigs")
            row = [
                rec.accession,
                rec.bioproject or st.bioproject,
                st.status,
                str(st.http) if st.http else "",
                str(st.n_contigs) if scanned else "",
                _fmt(st.total_bp / 1e6, 3) if scanned else "",
                _fmt(st.mapped_pct, 2) if scanned else "",
                _fmt(st.yield_pct, 2) if scanned else "",
                _fmt(st.divergence, 4) if scanned and st.divergence == st.divergence else "",
                str(st.n_tiles) if scanned else "",
                _fmt(st.tile_mass, 1) if scanned else "",
                str(st.n_spliced) if scanned else "",
                str(len(st.introns)) if scanned else "",
                str(st.rank) if st.rank else "",
                _fmt(st.gain, 3) if st.rank else "",
                _fmt(cum_s.get(rec.accession), 3) if st.rank else "",
            ]
            f.write("\t".join(row) + "\n")


def _samtools_merge(chunk_bams: List[Path], out: Path, threads: int, samtools: str = "samtools") -> None:
    if not chunk_bams:
        return
    if len(chunk_bams) == 1:
        shutil.copyfile(chunk_bams[0], out)
    else:
        cmd = [samtools, "merge", "-f", "-@", str(max(1, threads)), "-o", str(out),
               *[str(b) for b in chunk_bams]]
        log.info("samtools merge -> %s", out)
        subprocess.run(cmd, check=True)
    # BAI where possible (widest viewer support); CSI for chromosomes over
    # 512 Mbp, which BAI cannot hold (wheat 3B)
    try:
        subprocess.run([samtools, "index", str(out)], check=True, capture_output=True)
    except (subprocess.CalledProcessError, OSError):
        try:
            subprocess.run([samtools, "index", "-c", str(out)], check=True)
        except (subprocess.CalledProcessError, OSError) as e:  # pragma: no cover
            log.warning("samtools index %s failed: %s", out, e)


def run_logan(cfg: LoganConfig) -> int:
    """Run the whole Logan pre-screen. Returns a process exit code.

    0 = success, 3 = no run passed the quality gate, 4 = Logan S3 unreachable
    (every HEAD request errored).
    """
    t_start = time.monotonic()
    timings: Dict[str, float] = {}
    ldir = cfg.logan_dir
    contigs_dir = ldir / "contigs"
    runs_dir = ldir / "runs"
    bams_dir = ldir / "bams" if cfg.write_bam else ldir / "tmp"
    check_sequence_lengths(cfg.genome)
    for d in (ldir, contigs_dir, runs_dir, bams_dir):
        d.mkdir(parents=True, exist_ok=True)
    random.seed(cfg.seed)

    # ---- candidates -----------------------------------------------------
    records = read_runlist(cfg.runlist)
    rec_by_acc = {r.accession: r for r in records}
    order = sample_candidates(records, 0)  # full deterministic order
    n_target = cfg.max_candidates if cfg.max_candidates > 0 else len(order)
    log.info("%d runs in runlist; sampling up to %d candidates from %d bioprojects",
             len(records), n_target, len({r.bioproject for r in records if r.bioproject}))

    # ---- availability with top-up ----------------------------------------
    t0 = time.monotonic()
    avail: Dict[str, Tuple[str, int]] = {}
    cache_path = ldir / "availability.tsv"
    pos = 0
    n_available = 0
    while n_available < n_target and pos < len(order):
        batch = [r.accession for r in order[pos:pos + max(n_target - n_available, 1)]]
        pos += len(batch)
        res = check_availability(batch, cache_path, workers=cfg.download_workers,
                                 timeout=cfg.timeout, retries=cfg.retries)
        avail.update(res)
        n_available = sum(1 for a in avail.values() if a[0] == "available")
    timings["head"] = time.monotonic() - t0
    n_absent = sum(1 for a in avail.values() if a[0] == "absent")
    n_error = sum(1 for a in avail.values() if a[0] == "error")
    log.info("Availability: %d checked, %d available, %d absent, %d error (%.1f s)",
             len(avail), n_available, n_absent, n_error, timings["head"])
    if avail and n_error == len(avail):
        log.error("Logan S3 unreachable: all %d HEAD requests failed", n_error)
        return 4

    stats: Dict[str, LoganRunStats] = {}
    for acc, (st, code) in avail.items():
        if st != "available":
            stats[acc] = LoganRunStats(acc=acc, status=st, http=code,
                                       bioproject=rec_by_acc[acc].bioproject)
    todo = [a for a in (r.accession for r in order) if avail.get(a, ("", 0))[0] == "available"]

    # ---- resume: skip runs with saved stats ------------------------------
    pending: List[str] = []
    for acc in todo:
        st = load_run_stats(acc, runs_dir)
        if st is not None:
            st.bioproject = rec_by_acc[acc].bioproject
            st.status = "pending"
            st.rank, st.gain = 0, 0.0
            stats[acc] = st
        else:
            pending.append(acc)
    if len(todo) - len(pending):
        log.info("Resume: %d runs already scanned, %d to process", len(todo) - len(pending), len(pending))

    # ---- index ---------------------------------------------------------
    mmi = cfg.mmi
    if pending:
        if mmi is None:
            cand = ldir / "genome" / "mm2idx.mmi"
            if cand.is_file():
                log.info("Reusing minimap2 index %s", cand)
                mmi = cand
            else:
                t0 = time.monotonic()
                mmi = build_minimap2_index(cfg.genome, ldir / "genome", cfg.threads,
                                           part_bases=single_part_bases(cfg.genome))
                timings["index"] = time.monotonic() - t0
        elif not Path(mmi).is_file():
            raise FileNotFoundError(f"minimap2 index not found: {mmi}")

    # ---- download / align / scan pipeline ---------------------------------
    dl_bytes = 0
    dl_seconds = 0.0
    timings["align"] = 0.0
    timings["scan"] = 0.0
    chunk_bams: List[Path] = []
    chunk_no = 0
    t_pipe = time.monotonic()

    # The scan is single-threaded Python and takes about half as long as the
    # alignment of a chunk. By default minimap2's SAM is piped straight into
    # a scanner process (no BAM, no sort, nothing on disk), so a chunk is
    # scanned while it aligns and only a short tail remains when minimap2
    # exits. With --logan-bam the chunk is written as a sorted BAM first and
    # scanned from disk. Results are collected in submission order.
    scan_ex: Optional[cf.ProcessPoolExecutor] = None
    mp_ctx = mp.get_context("spawn")
    groups = 1
    if cfg.scan_workers > 0 and not cfg.write_bam:
        if cfg.align_groups is None:
            cfg.align_groups = auto_align_groups(cfg.threads)
        groups = max(1, int(cfg.align_groups))
        if groups > 1 and mmi is not None:
            fit = max_groups_for_memory(Path(mmi))
            if groups > fit:
                log.info("minimap2 processes per chunk %d -> %d: each loads its own "
                         "copy of the index (%s)", groups, fit, mmi)
                groups = fit
    if pending and mmi is not None:
        warn_if_index_does_not_fit(Path(mmi))
    n_scanners = max(cfg.scan_workers, groups) if cfg.scan_workers > 0 else 0
    if n_scanners > 0:
        scan_ex = cf.ProcessPoolExecutor(max_workers=n_scanners, mp_context=mp_ctx)
    # (future, contigs of the chunk, BAM or minimap2 log path, pipe read end)
    scan_pending: List[Tuple[cf.Future, List[LoganContigs], Path, object]] = []
    timings["scan_wait"] = 0.0  # time the main thread blocked on scan results

    # Thread budget: minimap2 runs beside the scanner processes (each scans
    # a chunk with one core, streamed or from disk) and beside the main
    # thread plus the download threads (network-bound; zstd decompression at
    # ~8 MB/s is negligible), which share one core. The index build and the
    # final LOGAN.bam merge run alone and keep all threads.
    mm2_threads = reserve_threads(cfg.threads, max(1, n_scanners) + 1)
    group_threads = max(1, mm2_threads // groups)
    log.info("Thread budget (--threads %d): minimap2 %d (%d x %d per chunk), scanners %d, "
             "main + downloads 1", cfg.threads, group_threads * groups, groups, group_threads,
             max(1, n_scanners))

    def collect_scan() -> None:
        fut, ready, path, conn = scan_pending.pop(0)
        t0 = time.monotonic()
        res, t_sc = fut.result()
        timings["scan_wait"] += time.monotonic() - t0
        timings["scan"] += t_sc
        if conn is not None and not conn.closed:
            conn.close()
        for acc, st in res.items():
            st.bioproject = rec_by_acc[acc].bioproject
            save_run_stats(st, runs_dir)
            stats[acc] = st
            log.info("  %-14s contigs=%6d mapped=%5.1f%% yield=%5.1f%% tiles=%5d spliced=%5d introns=%5d",
                     acc, st.n_contigs, st.mapped_pct, st.yield_pct, st.n_tiles, st.n_spliced,
                     len(st.introns))
        if cfg.write_bam:
            chunk_bams.append(path)
        else:
            path.unlink(missing_ok=True)  # the minimap2 log of a streamed chunk
        if not cfg.keep_contigs:
            for c in ready:
                c.fasta.unlink(missing_ok=True)
                (c.fasta.parent / f"{c.acc}.ka.npy").unlink(missing_ok=True)

    mem_groups = [groups]   # group count of the last chunk, to log changes only

    def process_chunk(ready: List[LoganContigs]) -> None:
        nonlocal chunk_no
        chunk_no += 1
        fastas = [c.fasta for c in ready]
        args = ({c.acc: c.ka for c in ready}, {c.acc: c.n_contigs for c in ready},
                {c.acc: c.total_bp for c in ready})
        kwargs = dict(tile_size=cfg.tile_size, ka_cap=cfg.ka_cap)
        t0 = time.monotonic()
        if cfg.write_bam:
            bam = bams_dir / f"chunk_{chunk_no:04d}.bam"
            align_contigs_minimap2(fastas, index=mmi, out_bam=bam,
                                   threads=mm2_threads, max_intron=cfg.max_intron,
                                   sort_threads=max(1, min(4, mm2_threads - 1)))
            t_al = time.monotonic() - t0
            if scan_ex is None:
                fut: cf.Future = cf.Future()
                fut.set_result(_scan_chunk_timed(bam, *args, **kwargs))
            else:
                fut = scan_ex.submit(_scan_chunk_timed, bam, *args, **kwargs)
            scan_pending.append((fut, ready, bam, None))
            summary = _minimap2_summary(bam.with_suffix(".minimap2.err"))
        else:
            # Memory is checked again before every chunk (the previous chunk's
            # minimap2 processes have exited): other jobs on the node may have
            # grown since the start. Fewer groups give the same results.
            n_groups = groups
            if groups > 1:
                n_groups = min(groups, max_groups_for_memory(Path(mmi)))
                if n_groups != mem_groups[0]:
                    log.info("chunk %d: %d minimap2 process(es) (%d configured), as many "
                             "copies of the index as fit in memory", chunk_no, n_groups, groups)
                    mem_groups[0] = n_groups
            # Balance the runs over the groups by contig bases (largest first).
            parts: List[List[LoganContigs]] = [[] for _ in range(min(n_groups, len(ready)))]
            load = [0] * len(parts)
            for c in sorted(ready, key=lambda c: c.total_bp, reverse=True):
                i = load.index(min(load))
                parts[i].append(c)
                load[i] += c.total_bp
            launched = []   # (proc, future, part, err, read end)
            try:
                for gi, part in enumerate(parts):
                    suffix = f"_g{gi}" if len(parts) > 1 else ""
                    err = bams_dir / f"chunk_{chunk_no:04d}{suffix}.minimap2.err"
                    pargs = ({c.acc: c.ka for c in part}, {c.acc: c.n_contigs for c in part},
                             {c.acc: c.total_bp for c in part})
                    r_conn, w_conn = mp_ctx.Pipe(duplex=False)
                    proc = start_contig_alignment([c.fasta for c in part], index=mmi,
                                                  stdout=w_conn.fileno(),
                                                  threads=max(1, mm2_threads // len(parts)),
                                                  max_intron=cfg.max_intron, log_path=err)
                    w_conn.close()  # minimap2 holds the only write end now
                    if scan_ex is None:
                        # inline scan (no pool): only ever one group here
                        fut = cf.Future()
                        fut.set_result(_scan_stream_timed(r_conn, *pargs, **kwargs))
                    else:
                        fut = scan_ex.submit(_scan_stream_timed, r_conn, *pargs, **kwargs)
                    launched.append((proc, fut, part, err, r_conn))
            except BaseException:
                for proc, *_ in launched:
                    proc.kill()
                    proc.wait()
                raise
            failed = None
            for proc, fut, part, err, r_conn in launched:
                rc = proc.wait()
                remove_split_files(err.with_suffix(".split"))  # split index only
                if rc != 0 and failed is None:
                    if fut.done() and fut.exception() is not None:
                        failed = fut.exception()  # the scanner died first; minimap2 got EPIPE
                    else:
                        failed = RuntimeError(f"minimap2 exited with status {rc}; see {err}")
            t_al = time.monotonic() - t0
            if failed is not None:
                raise failed
            for proc, fut, part, err, r_conn in launched:
                scan_pending.append((fut, part, err, r_conn))
            summary = "; ".join(x for x in (_minimap2_summary(e) for *_, e, _ in launched) if x)
        timings["align"] += t_al
        log.info("chunk %d: aligned %d runs in %.1f s%s", chunk_no, len(ready), t_al,
                 f" ({summary})" if summary else "")
        # bound the scans in flight, then reap what is done
        while len(scan_pending) > max(1, n_scanners) + groups or (
                scan_pending and scan_pending[0][0].done()):
            collect_scan()

    if pending:
        ready: List[LoganContigs] = []
        try:
            with cf.ThreadPoolExecutor(max_workers=max(1, cfg.download_workers)) as ex:
                futs = {ex.submit(download_contigs, acc, contigs_dir, timeout=cfg.timeout,
                                  retries=cfg.retries): acc for acc in pending}
                try:
                    for fut in cf.as_completed(futs):
                        acc = futs[fut]
                        try:
                            c = fut.result()
                        except LoganAbsentError as e:
                            log.warning("%s: absent at download time (HTTP %d)", acc, e.http)
                            stats[acc] = LoganRunStats(acc=acc, status="absent", http=e.http,
                                                       bioproject=rec_by_acc[acc].bioproject)
                            continue
                        except Exception as e:
                            log.error("%s: download failed: %s", acc, e)
                            stats[acc] = LoganRunStats(acc=acc, status="error",
                                                       http=int(getattr(e, "http", 0) or 0),
                                                       bioproject=rec_by_acc[acc].bioproject)
                            continue
                        dl_bytes += c.bytes_downloaded
                        dl_seconds += c.seconds
                        ready.append(c)
                        if len(ready) >= cfg.chunk_runs:
                            process_chunk(ready)
                            ready = []
                except BaseException:
                    ex.shutdown(wait=False, cancel_futures=True)
                    raise
            if ready:
                process_chunk(ready)
            while scan_pending:
                collect_scan()
        except BaseException:
            if scan_ex is not None:
                scan_ex.shutdown(wait=False, cancel_futures=True)
            raise
    if scan_ex is not None:
        scan_ex.shutdown(wait=True)
    timings["pipeline"] = time.monotonic() - t_pipe
    timings["download_thread_seconds"] = dl_seconds

    if cfg.write_bam and chunk_bams:
        t0 = time.monotonic()
        _samtools_merge(chunk_bams, ldir / "LOGAN.bam", cfg.threads)
        timings["merge"] = time.monotonic() - t0
    elif not cfg.write_bam:
        shutil.rmtree(bams_dir, ignore_errors=True)

    # ---- gate + selection --------------------------------------------------
    apply_gate(stats, min_contigs=cfg.min_contigs, min_tiles_frac=cfg.min_tiles_frac,
               max_divergence=cfg.max_divergence)
    for acc in rec_by_acc:
        if acc not in stats:
            stats[acc] = LoganRunStats(acc=acc, status="unsampled",
                                       bioproject=rec_by_acc[acc].bioproject)
    for st in stats.values():
        if st.status in ("accepted", "rejected", "too_few_contigs"):
            save_run_stats(st, runs_dir)

    from varus.select import greedy_select
    t0 = time.monotonic()
    scale = cfg.prior_batches * cfg.batch_size * 0.5
    curve = greedy_select(stats, top_k=cfg.select_top, scale=scale)
    timings["select"] = time.monotonic() - t0
    for acc, _, _, _ in curve:
        save_run_stats(stats[acc], runs_dir)

    # ---- outputs -------------------------------------------------------------
    _write_ranking(ldir / "LoganRanking.tsv", order, stats, curve)

    accepted = [s for s in stats.values() if s.status == "accepted"]
    ranked = sorted((s for s in accepted if s.rank), key=lambda s: s.rank)
    unranked = sorted((s for s in accepted if not s.rank), key=lambda s: (-s.n_tiles, s.acc))
    never = [r for r in records if stats[r.accession].status in ("absent", "error", "unsampled", "too_few_contigs")]
    out_records = [rec_by_acc[s.acc] for s in ranked] + [rec_by_acc[s.acc] for s in unranked] + never
    n_written = write_runlist(out_records, cfg.outdir / "Runlist.logan.tsv")
    log.info("Wrote %d runs to %s (%d ranked, %d accepted-unranked, %d unprocessed)",
             n_written, cfg.outdir / "Runlist.logan.tsv", len(ranked), len(unranked), len(never))

    # introns of accepted runs -> strand -> gff / splice sites / junc bed
    union: Dict[IntronKey, float] = {}
    for s in accepted:
        for k, w in s.introns.items():
            union[k] = union.get(k, 0.0) + w
    raw = IntronCounts({k: max(1, int(round(w))) for k, w in union.items()})
    assigner = StrandAssigner(cfg.genome)
    try:
        stranded, _ = assigner.assign_new(raw)
    finally:
        assigner.close()
    n_gff = write_introns_gff(stranded, ldir / "logan_introns.gff")
    write_hisat2_splice_sites(stranded, ldir / "logan.splice_sites")
    write_minimap2_junc_bed(stranded, ldir / "logan.junc.bed")
    log.info("Introns: %d raw from %d accepted runs, %d stranded -> %s",
             len(raw), len(accepted), n_gff, ldir / "logan_introns.gff")

    with gzip.open(ldir / "logan_tiles.tsv.gz", "wt", encoding="utf-8") as f:
        f.write("acc\tchrom\ttile_idx\tweight\n")
        for s in sorted(accepted, key=lambda s: (s.rank or 10**9, s.acc)):
            for (chrom, idx), w in sorted(s.tiles.items()):
                f.write(f"{s.acc}\t{chrom}\t{idx}\t{w:.3f}\n")

    counts: Dict[str, int] = {}
    for s in stats.values():
        counts[s.status] = counts.get(s.status, 0) + 1
    timings["total"] = time.monotonic() - t_start
    params = {k: (str(v) if isinstance(v, Path) else v) for k, v in asdict(cfg).items()}
    summary = {
        "params": params,
        "n_runlist": len(records),
        "n_candidates_checked": len(avail),
        "counts": counts,
        "timings_s": {k: round(v, 3) for k, v in timings.items()},
        "download_mb": round(dl_bytes / 1e6, 3),
        "download_mb_per_s": round(dl_bytes / 1e6 / dl_seconds, 3) if dl_seconds > 0 else 0.0,
        "scale": scale,
        "n_introns_stranded": n_gff,
        "coverage_curve": [
            {"rank": i + 1, "acc": acc, "gain": round(g, 4), "cumulative_S": round(S, 4),
             "cumulative_tiles": nt}
            for i, (acc, g, S, nt) in enumerate(curve)
        ],
    }
    (ldir / "logan_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    log.info("Status counts: %s; total %.1f s", counts, timings["total"])
    if not accepted:
        log.error("No run passed the quality gate")
        return 3
    return 0

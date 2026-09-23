"""``varus logan``: pre-screen SRA runs with Logan contig assemblies.

Logan (Chikhi et al. 2024) provides per-run assemblies of essentially all
public SRA runs up to the end of 2023. Instead of blindly downloading reads
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

from varus.align import align_contigs_minimap2
from varus.index import build_minimap2_index
from varus.introns import IntronCounts, IntronKey, iter_introns, write_introns_gff
from varus.runlist import RunRecord, write_runlist
from varus.strand import StrandAssigner, write_hisat2_splice_sites, write_minimap2_junc_bed
from varus.tiles import Tile

log = logging.getLogger(__name__)

LOGAN_CONTIGS_URL = "https://s3.amazonaws.com/logan-pub/c/{acc}/{acc}.contigs.fa.zst"
LOGAN_UNITIGS_URL = "https://s3.amazonaws.com/logan-pub/u/{acc}/{acc}.unitigs.fa.zst"

TILE_WEIGHTS = ("unit", "ka", "ka_len")
STATUSES = ("pending", "accepted", "rejected", "too_few_contigs", "absent", "error", "unsampled")

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
    download_workers: int = 8
    max_candidates: int = 500
    chunk_runs: int = 10
    max_intron: int = 20_000
    min_contigs: int = 100
    min_tiles_frac: float = 0.10
    ka_cap: float = 50.0
    tile_weight: str = "ka_len"  # unit | ka | ka_len
    tile_size: int = 5000
    select_top: int = 50
    batch_size: int = 50_000
    prior_batches: float = 1.0
    keep_contigs: bool = False
    write_bam: bool = False
    seed: int = 1
    timeout: float = 60.0
    retries: int = 4
    # Recorded in logan_summary.json for the downstream ``varus run`` stage;
    # contig alignment itself is platform-agnostic.
    longreads: bool = False

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
    bam: Path,
    ka_by_acc: Dict[str, np.ndarray],
    n_contigs_by_acc: Dict[str, int],
    total_bp_by_acc: Dict[str, int],
    *,
    tile_size: int,
    ka_cap: float,
    tile_weight: str,
) -> Dict[str, LoganRunStats]:
    """One pysam pass over a chunk BAM; per-run tile and intron weights.

    Only primary alignments are counted (unmapped, secondary and
    supplementary records are skipped). The run is recovered from the query
    name ``<ACC>_<i>``. Tile weight per contig:

    * ``unit``   - 1
    * ``ka``     - ``min(ka, ka_cap)`` (NaN -> 1)
    * ``ka_len`` - ``min(ka, ka_cap) * max(1, aligned_ref_len / 150)``

    Introns (CIGAR ``N``) are weighted with ``min(ka, ka_cap)`` regardless of
    ``tile_weight``. Every accession in ``ka_by_acc`` gets a stats entry even
    if none of its contigs aligned.
    """
    if tile_weight not in TILE_WEIGHTS:
        raise ValueError(f"tile_weight must be one of {TILE_WEIGHTS}, got {tile_weight!r}")
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
    with pysam.AlignmentFile(str(bam), "rb") as fh:
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
            if tile_weight == "unit":
                w = 1.0
            elif tile_weight == "ka":
                w = wk
            else:
                ref_len = read.reference_end - read.reference_start if read.reference_end is not None else 0
                w = wk * max(1.0, ref_len / 150.0)
            tile: Tile = (read.reference_name, (read.reference_start + 1) // tile_size)
            st.tiles[tile] = st.tiles.get(tile, 0.0) + w
            st.n_aligned += 1
            st.aligned_bp += int(read.query_alignment_length)
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
    return stats


def apply_gate(
    stats: Dict[str, LoganRunStats],
    *,
    min_contigs: int,
    min_tiles_frac: float,
) -> None:
    """Set ``status`` of every scanned run in place.

    * ``too_few_contigs`` if ``n_contigs < min_contigs``,
    * otherwise ``accepted`` if ``n_tiles >= min_tiles_frac * max_tiles``
      (``max_tiles`` = best ``n_tiles`` among runs with enough contigs) and
      ``n_tiles > 0``, else ``rejected``.

    Runs whose status is ``absent`` / ``error`` / ``unsampled`` are untouched.
    """
    scannable = [s for s in stats.values()
                 if s.status in ("pending", "accepted", "rejected", "too_few_contigs")]
    eligible = [s for s in scannable if s.n_contigs >= min_contigs]
    max_tiles = max((s.n_tiles for s in eligible), default=0)
    thr = min_tiles_frac * max_tiles
    n_acc = n_rej = n_few = 0
    for s in scannable:
        if s.n_contigs < min_contigs:
            s.status = "too_few_contigs"
            n_few += 1
        elif s.n_tiles > 0 and s.n_tiles >= thr:
            s.status = "accepted"
            n_acc += 1
        else:
            s.status = "rejected"
            n_rej += 1
    log.info("Gate: max_tiles=%d threshold=%.1f -> %d accepted, %d rejected, %d too few contigs",
             max_tiles, thr, n_acc, n_rej, n_few)


# ---------------------------------------------------------------------------
# Persistence of per-run stats (resume)
# ---------------------------------------------------------------------------

_SCALAR_FIELDS = ("acc", "n_contigs", "total_bp", "n_aligned", "aligned_bp", "mapped_pct",
                  "mass_total", "mass_mapped", "yield_pct",
                  "n_spliced", "n_tiles", "tile_mass", "status", "rank", "gain", "http",
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
    scalars = {k: getattr(stats, k) for k in _SCALAR_FIELDS}
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
        if k in scalars:
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
            "yield_pct", "n_tiles", "tile_mass", "n_spliced", "n_introns", "selected_rank",
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
                f.write(f"{rec.accession}\t{rec.bioproject}\tunsampled\t" + "\t" * 12 + "\n")
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
    try:
        subprocess.run([samtools, "index", str(out)], check=True)
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
    bams_dir = ldir / "bams" if cfg.write_bam else ldir / "tmp_bams"
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
                mmi = build_minimap2_index(cfg.genome, ldir / "genome", cfg.threads)
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

    def process_chunk(ready: List[LoganContigs]) -> None:
        nonlocal chunk_no
        chunk_no += 1
        bam = bams_dir / f"chunk_{chunk_no:04d}.bam"
        t0 = time.monotonic()
        align_contigs_minimap2([c.fasta for c in ready], index=mmi, out_bam=bam,
                               threads=cfg.threads, max_intron=cfg.max_intron)
        t_al = time.monotonic() - t0
        timings["align"] += t_al
        log.info("chunk %d: aligned %d runs in %.1f s", chunk_no, len(ready), t_al)
        t0 = time.monotonic()
        res = scan_chunk_bam(
            bam,
            {c.acc: c.ka for c in ready},
            {c.acc: c.n_contigs for c in ready},
            {c.acc: c.total_bp for c in ready},
            tile_size=cfg.tile_size, ka_cap=cfg.ka_cap, tile_weight=cfg.tile_weight,
        )
        for acc, st in res.items():
            st.bioproject = rec_by_acc[acc].bioproject
            save_run_stats(st, runs_dir)
            stats[acc] = st
            log.info("  %-14s contigs=%6d mapped=%5.1f%% yield=%5.1f%% tiles=%5d spliced=%5d introns=%5d",
                     acc, st.n_contigs, st.mapped_pct, st.yield_pct, st.n_tiles, st.n_spliced,
                     len(st.introns))
        timings["scan"] += time.monotonic() - t0
        if cfg.write_bam:
            chunk_bams.append(bam)
        else:
            bam.unlink(missing_ok=True)
            err = bam.with_suffix(".minimap2.err")
            err.unlink(missing_ok=True)
        if not cfg.keep_contigs:
            for c in ready:
                c.fasta.unlink(missing_ok=True)
                (c.fasta.parent / f"{c.acc}.ka.npy").unlink(missing_ok=True)

    if pending:
        ready: List[LoganContigs] = []
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
    timings["pipeline"] = time.monotonic() - t_pipe
    timings["download_thread_seconds"] = dl_seconds

    if cfg.write_bam and chunk_bams:
        t0 = time.monotonic()
        _samtools_merge(chunk_bams, ldir / "LOGAN.bam", cfg.threads)
        timings["merge"] = time.monotonic() - t0
    elif not cfg.write_bam:
        shutil.rmtree(bams_dir, ignore_errors=True)

    # ---- gate + selection --------------------------------------------------
    apply_gate(stats, min_contigs=cfg.min_contigs, min_tiles_frac=cfg.min_tiles_frac)
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

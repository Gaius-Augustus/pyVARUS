"""Alignment of one batch (HISAT2 for short reads, minimap2 for long reads).

The HISAT2 path improves on ``legacy/Implementation/src/HISAT_Aligner.cpp``:

* No SAM intermediate. HISAT2's stdout is piped through ``samtools sort -O BAM``
  directly, eliminating one full read+write of the alignments.
* No separate ``samtools sort`` step afterwards.
* ``Log.final.out`` (HISAT2's stderr summary) is captured to a file in the
  batch directory so the legacy quality-parsing logic still works on it.

The minimap2 path follows the same pipe-to-``samtools sort`` pattern; minimap2
emits no end-of-run summary, so the quality gate is computed by scanning the
sorted BAM with pysam (see :func:`count_minimap2_quality`).

The output is ``<batch_dir>/Aligned.out.bam`` (coordinate-sorted), matching the
filename the legacy controller already expects after its own post-conversion.
"""

from __future__ import annotations

import functools
import logging
import shutil
import struct
import subprocess
from dataclasses import dataclass
from pathlib import Path

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class AlignmentResult:
    bam: Path
    log: Path


def _require(tool: str) -> None:
    if shutil.which(tool) is None:
        raise RuntimeError(f"{tool} not found on PATH")


def reserve_threads(total: int, reserved: int) -> int:
    """Threads left for the main tool after ``reserved`` for concurrent work.

    Reservations are capped at a quarter of ``total``: with very few threads
    (e.g. 2) halving the aligner costs more than briefly oversubscribing.
    """
    total = max(1, int(total))
    return max(1, total - min(max(0, reserved), total // 4))


def align_batch_hisat2(
    r1: Path,
    r2: Path | None,
    *,
    index_prefix: Path,
    batch_dir: Path,
    threads: int = 4,
    intron_db: Path | None = None,
    hisat2: str = "hisat2",
    samtools: str = "samtools",
    mm: bool = True,
    keep_unaligned: bool = False,
    sort_compression: int | None = 1,
    sort_threads: int | None = None,
) -> AlignmentResult:
    """Align one batch with HISAT2; write a sorted BAM.

    Parameters
    ----------
    r1, r2
        FASTA files. ``r2`` is ``None`` for single-end reads.
    index_prefix
        HISAT2 index prefix (e.g. ``<outdir>/genome/hisatidx``).
    batch_dir
        Output directory; ``Aligned.out.bam`` and ``Log.final.out`` go there.
    intron_db
        Optional path to a known-splice-site file in HISAT2's tab format
        (``--known-splicesite-infile``). If absent, alignment proceeds without
        a splice DB, like the very first batch in the legacy code.
    mm
        Pass ``--mm`` so HISAT2 memory-maps the index. Successive (and
        concurrent) invocations then share the OS page cache instead of each
        re-reading the index into private memory. Disable on file systems
        where mmap is slow (some NFS setups).
    keep_unaligned
        By default ``--no-unal`` drops unaligned reads from the BAM; nothing
        downstream (tile counting, intron extraction, BRAKER) uses them and
        they inflate every per-batch BAM, the sort, and the final merge.
        The quality gate is parsed from HISAT2's log, so it is unaffected.
    sort_compression
        ``samtools sort -l`` level for the per-batch BAM. Per-batch BAMs are
        re-compressed by the final merge, so a low level is cheapest overall.
        ``None`` keeps the samtools default.
    sort_threads
        ``samtools sort -@``; default ``threads - 1``.
    """
    _require(hisat2)
    _require(samtools)
    batch_dir.mkdir(parents=True, exist_ok=True)
    bam_out = batch_dir / "Aligned.out.bam"
    log_out = batch_dir / "Log.final.out"

    hisat_cmd: list[str] = [
        hisat2,
        "-p", str(threads),
        "-f",
        "-x", str(index_prefix),
    ]
    if mm:
        hisat_cmd.append("--mm")
    if not keep_unaligned:
        hisat_cmd.append("--no-unal")
    if r2 is None:
        hisat_cmd += ["-U", str(r1)]
    else:
        hisat_cmd += ["-1", str(r1), "-2", str(r2)]

    if intron_db is not None and intron_db.is_file():
        hisat_cmd += ["--known-splicesite-infile", str(intron_db)]

    sort_cmd = [
        samtools, "sort",
        "-@", str(sort_threads if sort_threads else max(1, threads - 1)),
        "-O", "BAM",
    ]
    if sort_compression is not None:
        sort_cmd += ["-l", str(int(sort_compression))]
    sort_cmd += ["-o", str(bam_out)]

    log.info("HISAT2 | samtools sort -> %s", bam_out)
    log.debug("hisat2 cmd: %s", " ".join(hisat_cmd))
    log.debug("samtools cmd: %s", " ".join(sort_cmd))

    with log_out.open("wb") as logf:
        hisat_proc = subprocess.Popen(
            hisat_cmd, stdout=subprocess.PIPE, stderr=logf
        )
        try:
            sort_proc = subprocess.Popen(
                sort_cmd, stdin=hisat_proc.stdout, stdout=subprocess.DEVNULL
            )
            # Allow hisat2 to receive SIGPIPE if samtools dies first.
            assert hisat_proc.stdout is not None
            hisat_proc.stdout.close()
            sort_rc = sort_proc.wait()
        finally:
            hisat_rc = hisat_proc.wait()

    if hisat_rc != 0:
        raise RuntimeError(
            f"hisat2 exited with status {hisat_rc}; see {log_out}"
        )
    if sort_rc != 0:
        raise RuntimeError(f"samtools sort exited with status {sort_rc}")

    return AlignmentResult(bam=bam_out, log=log_out)


LONGREAD_PRESETS: dict[str, list[str]] = {
    "pacbio": ["-ax", "splice"],
    "ont": ["-ax", "splice", "-uf", "-k14"],
}

# Map SRA platform identifiers (as parsed from esummary <Instrument>) to the
# minimap2 preset key used by :func:`align_batch_minimap2`.
_PLATFORM_TO_PRESET: dict[str, str] = {
    "PACBIO_SMRT": "pacbio",
    "OXFORD_NANOPORE": "ont",
}


def preset_for_platform(platform: str) -> str:
    """Return the minimap2 preset for an SRA platform string.

    Falls back to ``'pacbio'`` (the dominant long-read RNA-seq submission) and
    logs a warning when the platform is missing or unrecognised — the user can
    correct this by editing the ``platform`` column in ``Runlist.tsv``.
    """
    key = (platform or "").upper().strip()
    preset = _PLATFORM_TO_PRESET.get(key)
    if preset is not None:
        return preset
    log.warning(
        "Unknown / missing platform %r; defaulting to minimap2 'pacbio' preset.",
        platform,
    )
    return "pacbio"


def align_batch_minimap2(
    reads: Path,
    *,
    index: Path,
    batch_dir: Path,
    threads: int = 4,
    preset: str = "pacbio",
    junc_bed: Path | None = None,
    minimap2: str = "minimap2",
    samtools: str = "samtools",
    sort_threads: int | None = None,
) -> AlignmentResult:
    """Align one batch of long reads with minimap2; write a sorted BAM.

    Parameters
    ----------
    reads
        FASTA file. Long-read SRA runs are single-end; ``r2`` is never used.
    index
        ``.mmi`` index file from :func:`varus.index.build_minimap2_index`, or
        a genome FASTA (minimap2 will index on the fly).
    batch_dir
        Output directory; ``Aligned.out.bam`` and ``Log.minimap2.err`` go there.
    preset
        ``'pacbio'`` (Iso-Seq / HiFi, ``-ax splice``) or ``'ont'``
        (direct-RNA Nanopore, ``-ax splice -uf -k14``).
    junc_bed
        Optional BED12 of known junctions for minimap2 ``--junc-bed``. If
        absent, alignment proceeds without the hint, like the first batch.
    """
    if preset not in LONGREAD_PRESETS:
        raise ValueError(
            f"unknown long-read preset {preset!r}; "
            f"expected one of {list(LONGREAD_PRESETS)}"
        )
    _require(minimap2)
    _require(samtools)
    batch_dir.mkdir(parents=True, exist_ok=True)
    bam_out = batch_dir / "Aligned.out.bam"
    log_out = batch_dir / "Log.minimap2.err"

    mm2_cmd: list[str] = [
        minimap2,
        "-t", str(threads),
        *LONGREAD_PRESETS[preset],
    ]
    if junc_bed is not None and junc_bed.is_file():
        mm2_cmd += ["--junc-bed", str(junc_bed)]
    split = minimap2_split_prefix(index, batch_dir / "minimap2.split")
    if split is not None:
        mm2_cmd += ["--split-prefix", str(split)]
    mm2_cmd += [str(index), str(reads)]

    sort_cmd = [
        samtools, "sort",
        "-@", str(sort_threads if sort_threads else max(1, threads - 1)),
        "-O", "BAM",
        "-o", str(bam_out),
    ]

    log.info("minimap2 (%s) | samtools sort -> %s", preset, bam_out)
    log.debug("minimap2 cmd: %s", " ".join(mm2_cmd))
    log.debug("samtools cmd: %s", " ".join(sort_cmd))

    with log_out.open("wb") as logf:
        mm2_proc = subprocess.Popen(
            mm2_cmd, stdout=subprocess.PIPE, stderr=logf
        )
        try:
            sort_proc = subprocess.Popen(
                sort_cmd, stdin=mm2_proc.stdout, stdout=subprocess.DEVNULL
            )
            assert mm2_proc.stdout is not None
            mm2_proc.stdout.close()
            sort_rc = sort_proc.wait()
        finally:
            mm2_rc = mm2_proc.wait()

    if mm2_rc != 0:
        raise RuntimeError(
            f"minimap2 exited with status {mm2_rc}; see {log_out}"
        )
    if sort_rc != 0:
        raise RuntimeError(f"samtools sort exited with status {sort_rc}")

    return AlignmentResult(bam=bam_out, log=log_out)


_MMI_MAGIC = b"MMI\x02"
_MM_I_NO_SEQ = 0x2
# A genome FASTA passed as minimap2 index is indexed on the fly and split
# every 8 Gbp (-I); above this file size (gzip ~4x) it may be split.
_FASTA_MAY_SPLIT_BYTES = 1_500_000_000


def minimap2_index_parts(index: Path) -> list[int] | None:
    """Byte size of every part of a minimap2 ``.mmi`` index.

    minimap2 splits an index every ~8 Gbp of genome (``-I 8G``); the parts
    are complete indexes dumped one after the other (``mm_idx_dump`` in
    minimap2's index.c: magic, w/k/b/n_seq/flag, names and lengths, 2^b hash
    buckets, packed sequence). minimap2 holds one part in memory at a time.
    None if ``index`` is unreadable, not a ``.mmi`` or not laid out as expected.
    Cached per file: the walk touches pages all over the index (5 s for a
    cold 444 MB index on ceph, 0.1 s warm).
    """
    try:
        st = Path(index).stat()
    except OSError:
        return None
    parts = _index_parts(str(Path(index).resolve()), st.st_size, st.st_mtime)
    return list(parts) if parts is not None else None


@functools.lru_cache(maxsize=None)
def _index_parts(path: str, size: int, mtime: float) -> tuple[int, ...] | None:
    try:
        fh = open(path, "rb")
    except OSError:
        return None
    parts: list[int] = []
    u32 = struct.Struct("<I")
    with fh:
        start = 0
        while start < size:
            fh.seek(start)
            hdr = fh.read(24)
            if len(hdr) != 24 or hdr[:4] != _MMI_MAGIC:
                return None
            _w, _k, b, n_seq, flag = struct.unpack("<5I", hdr[4:])
            if b > 32:
                return None
            sum_len = 0
            for _ in range(n_seq):
                name_len = fh.read(1)
                if not name_len:
                    return None
                fh.seek(name_len[0], 1)
                x = fh.read(4)
                if len(x) != 4:
                    return None
                sum_len += u32.unpack(x)[0]
            for _ in range(1 << b):
                x = fh.read(4)
                if len(x) != 4:
                    return None
                fh.seek(8 * u32.unpack(x)[0], 1)      # b->p, uint64 each
                x = fh.read(4)
                if len(x) != 4:
                    return None
                fh.seek(16 * u32.unpack(x)[0], 1)     # hash (key, value) pairs
            if not flag & _MM_I_NO_SEQ:
                fh.seek(4 * ((sum_len + 7) // 8), 1)  # 4-bit packed sequence
            end = fh.tell()
            if end > size:
                return None
            parts.append(end - start)
            start = end
    return tuple(parts) or None


@functools.lru_cache(maxsize=None)
def _index_is_split(path: str, size: int, mtime: float) -> bool:
    index = Path(path)
    try:
        with index.open("rb") as fh:
            head = fh.read(4)
    except OSError:
        return False
    if head[:3] == _MMI_MAGIC[:3]:
        parts = minimap2_index_parts(index)
        if parts is None:
            log.warning("Could not read the part layout of the minimap2 index %s; running "
                        "minimap2 with --split-prefix, which is correct for any index but "
                        "writes its output only after all alignments are done.", index)
            return True
        if len(parts) > 1:
            log.info("minimap2 index %s has %d parts (genome > ~8 Gbp, largest part "
                     "%.1f GB in memory at a time); minimap2 runs with --split-prefix so "
                     "the parts' alignments are merged into one correct SAM.",
                     index, len(parts), max(parts) / 2**30)
            return True
        return False
    if size > _FASTA_MAY_SPLIT_BYTES:
        log.info("Genome FASTA %s used as minimap2 index is large enough to be split "
                 "(> 8 Gbp); minimap2 runs with --split-prefix.", index)
        return True
    return False


def minimap2_split_prefix(index: Path, prefix: Path) -> Path | None:
    """``prefix`` if ``index`` is (or may be) split into several parts, else None.

    Without ``--split-prefix`` minimap2 writes the alignments of every part
    separately: no ``@SQ`` header, every read once per part (unmapped in the
    parts without its locus) and a multi-mapper that spans two parts gets
    two primary alignments, each with MAPQ 60. With it, minimap2 keeps each
    part's alignments in ``<prefix>.NNNN.tmp`` files and merges them at the
    end; the output then matches a single-part index. It is only passed for
    split indexes because the merge holds all output back until the end.
    """
    try:
        st = Path(index).stat()
    except OSError:
        return None
    if not _index_is_split(str(Path(index).resolve()), st.st_size, st.st_mtime):
        return None
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    return Path(prefix)


def _contig_minimap2_cmd(
    queries: list[Path],
    index: Path,
    threads: int,
    max_intron: int,
    minimap2: str,
    mini_batch: str | None,
    split_prefix: Path | None = None,
) -> list[str]:
    cmd: list[str] = [
        minimap2,
        "-t", str(threads),
        "-ax", "splice",
        "--secondary=no",
        "-G", str(int(max_intron)),
    ]
    if mini_batch:
        cmd += ["-K", str(mini_batch)]
    if split_prefix is not None:
        cmd += ["--split-prefix", str(split_prefix)]
    cmd += [str(index), *[str(q) for q in queries]]
    return cmd


def start_contig_alignment(
    queries: list[Path],
    *,
    index: Path,
    stdout: int,
    threads: int = 4,
    max_intron: int = 20_000,
    log_path: Path,
    minimap2: str = "minimap2",
    mini_batch: str | None = None,
) -> subprocess.Popen:
    """Start ``minimap2 -ax splice`` on Logan contigs, SAM to file descriptor ``stdout``.

    The caller owns the descriptor (typically the write end of a pipe whose
    read end a scanner process consumes) and must close its own copy after
    this returns, so the reader sees EOF when minimap2 exits. ``mini_batch``
    sets minimap2's ``-K``; the default (None, minimap2's 500 Mb) is kept
    because ``-K 20M`` left 48-thread minimap2 at 10–25 busy cores and made
    the Drosophila stage 18 % slower (2026-09-24). The scanner runs in its
    own process, so streaming within a chunk buys nothing.

    Same alignment settings as :func:`align_contigs_minimap2`: ``-ax splice``
    without the ONT/PacBio error-model tweaks (contigs are consensus
    sequences), ``--secondary=no`` so every contig contributes at most one
    primary alignment per locus, ``-G`` capping the intron length. Several
    query FASTAs per call amortise the index load over a chunk of runs; read
    names keep their ``<ACC>_<i>`` prefix so the scanner can split by run.
    """
    _require(minimap2)
    if not queries:
        raise ValueError("start_contig_alignment: no query FASTA given")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    split = minimap2_split_prefix(index, log_path.with_suffix(".split"))
    cmd = _contig_minimap2_cmd(queries, index, threads, max_intron, minimap2, mini_batch, split)
    log.info("minimap2 splice (contigs, %d files) -> scanner", len(queries))
    log.debug("minimap2 cmd: %s", " ".join(cmd))
    with log_path.open("wb") as logf:
        return subprocess.Popen(cmd, stdout=stdout, stderr=logf)


def align_contigs_minimap2(
    queries: list[Path],
    *,
    index: Path,
    out_bam: Path,
    threads: int = 4,
    max_intron: int = 20_000,
    log_path: Path | None = None,
    minimap2: str = "minimap2",
    samtools: str = "samtools",
    sort_compression: int | None = 1,
    sort_threads: int | None = None,
) -> Path:
    """Spliced-align assembled contigs (Logan) to the genome; write a sorted BAM.

    Used only when the chunk BAMs are to be kept (``varus logan --logan-bam``);
    the default path streams SAM into the scanner via
    :func:`start_contig_alignment`. Settings as documented there.
    """
    _require(minimap2)
    _require(samtools)
    if not queries:
        raise ValueError("align_contigs_minimap2: no query FASTA given")
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    log_out = log_path or out_bam.with_suffix(".minimap2.err")

    split = minimap2_split_prefix(index, out_bam.with_suffix(".split"))
    mm2_cmd = _contig_minimap2_cmd(queries, index, threads, max_intron, minimap2, None, split)
    sort_cmd = [samtools, "sort", "-@",
                str(sort_threads if sort_threads else max(1, threads - 1)), "-O", "BAM"]
    if sort_compression is not None:
        sort_cmd += ["-l", str(int(sort_compression))]
    sort_cmd += ["-o", str(out_bam)]

    log.info("minimap2 splice (contigs, %d files) | samtools sort -> %s",
             len(queries), out_bam)
    log.debug("minimap2 cmd: %s", " ".join(mm2_cmd))
    with log_out.open("wb") as logf:
        mm2_proc = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=logf)
        try:
            sort_proc = subprocess.Popen(
                sort_cmd, stdin=mm2_proc.stdout, stdout=subprocess.DEVNULL
            )
            assert mm2_proc.stdout is not None
            mm2_proc.stdout.close()
            sort_rc = sort_proc.wait()
        finally:
            mm2_rc = mm2_proc.wait()
    if mm2_rc != 0:
        raise RuntimeError(f"minimap2 exited with status {mm2_rc}; see {log_out}")
    if sort_rc != 0:
        raise RuntimeError(f"samtools sort exited with status {sort_rc}")
    return out_bam


def count_minimap2_quality(
    bam_path: Path,
    *,
    min_mapq: int = 1,
) -> dict[str, float]:
    """Count primary, MAPQ ≥ ``min_mapq`` alignments in a minimap2 BAM.

    Returns the same dict shape as :func:`parse_hisat2_log`:
    ``{'num_uniq': float, 'uniq_pct': float}``.

    Note the denominator difference vs ``parse_hisat2_log``: HISAT2's parser
    divides by ``batch_size`` (input pairs/reads); here we divide by primary
    alignments observed. For long-read SRA spots-are-reads, this is the
    fraction of decoded reads that aligned uniquely — practically equivalent
    for the ``--min-uniq-pct`` quality gate but slightly differently calibrated.
    """
    if not bam_path.is_file():
        raise FileNotFoundError(bam_path)

    import pysam  # local import: pysam is in the optional [align] extra

    n_primary = 0
    n_uniq = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
            ):
                continue
            n_primary += 1
            if read.mapping_quality >= min_mapq:
                n_uniq += 1

    if n_primary == 0:
        return {"num_uniq": 0.0, "uniq_pct": 0.0}
    return {
        "num_uniq": float(n_uniq),
        "uniq_pct": 100.0 * n_uniq / n_primary,
    }


def parse_hisat2_log(log_path: Path, batch_size: int) -> dict[str, float]:
    """Parse HISAT2's ``Log.final.out`` for unique-alignment percentage.

    Mirrors ``HISAT_Aligner::checkQuality``: the first occurrence of either
    "aligned concordantly exactly 1 time" (paired) or "aligned exactly 1 time"
    (single-end) gives the unique alignment count.
    """
    if not log_path.is_file():
        raise FileNotFoundError(log_path)
    num_uniq: int | None = None
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.lstrip()
        if num_uniq is None and (
            "aligned concordantly exactly 1 time" in s
            or "aligned exactly 1 time" in s
        ):
            try:
                num_uniq = int(s.split()[0])
            except (ValueError, IndexError):
                continue
            break
    if num_uniq is None:
        return {"num_uniq": 0, "uniq_pct": 0.0}
    return {
        "num_uniq": float(num_uniq),
        "uniq_pct": 100.0 * num_uniq / batch_size,
    }

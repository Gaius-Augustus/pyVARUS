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

import logging
import shutil
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
        "-@", str(max(1, threads - 1)),
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
    mm2_cmd += [str(index), str(reads)]

    sort_cmd = [
        samtools, "sort",
        "-@", str(max(1, threads - 1)),
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


def _contig_minimap2_cmd(
    queries: list[Path],
    index: Path,
    threads: int,
    max_intron: int,
    minimap2: str,
    mini_batch: str | None,
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
    cmd = _contig_minimap2_cmd(queries, index, threads, max_intron, minimap2, mini_batch)
    log_path.parent.mkdir(parents=True, exist_ok=True)
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

    mm2_cmd = _contig_minimap2_cmd(queries, index, threads, max_intron, minimap2, None)
    sort_cmd = [samtools, "sort", "-@", str(max(1, threads - 1)), "-O", "BAM"]
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

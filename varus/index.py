"""Build a HISAT2 or minimap2 index for a genome FASTA.

Replaces the index-building branch of ``legacy/runVARUS.pl``.

For short reads we invoke ``hisat2-build`` directly. The output prefix matches
what the legacy aligner expected (``<outdir>/hisatidx``), so a ``varus run``
task can point at ``--index <outdir>`` and HISAT2 will resolve
``<outdir>/hisatidx.*.ht2``.

For long reads (``--longreads``) we invoke ``minimap2 -d`` with the ``splice``
preset. minimap2 can index on the fly, but pre-building a ``.mmi`` saves the
seed-table construction cost on every batch. The ``splice`` index works for
both PacBio Iso-Seq and ONT direct-RNA; only the alignment-time flags differ.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

from varus.align import minimap2_index_parts

log = logging.getLogger(__name__)

# BAM stores positions as int32, and minimap2 (2.31) cannot parse a longer
# FASTA record: it logs "failed to parse the first FASTA/FASTQ record" and
# leaves the sequence out of the index, so its reads silently go unmapped.
MAX_SEQUENCE_LENGTH = 2**31 - 1


def _fai_lengths(genome: Path) -> list[tuple[str, int]] | None:
    """(name, length) of every sequence, from ``<genome>.fai`` (built with
    pyfaidx if missing or older than the FASTA); None if that fails."""
    fai = Path(str(genome) + ".fai")
    try:
        if not fai.is_file() or fai.stat().st_mtime < Path(genome).stat().st_mtime:
            from pyfaidx import Faidx
            Faidx(str(genome)).close()
        out = []
        for line in fai.read_text().splitlines():
            name, length = line.split("\t")[:2]
            out.append((name, int(length)))
        return out
    except Exception as e:  # pyfaidx missing, gzip without bgzf, read-only dir
        log.warning("Could not index %s to check its sequence lengths (%s)", genome, e)
        return None


def check_sequence_lengths(genome: Path) -> None:
    """Raise ``ValueError`` if a sequence of ``genome`` is longer than
    2^31 - 1 bp, which neither BAM nor minimap2 can handle (axolotl and
    lungfish chromosomes are). Uses and, if needed, writes ``<genome>.fai``."""
    lengths = _fai_lengths(Path(genome))
    if not lengths:
        return
    too_long = [(n, ln) for n, ln in lengths if ln > MAX_SEQUENCE_LENGTH]
    if too_long:
        shown = ", ".join(f"{n} ({ln / 1e9:.2f} Gbp)" for n, ln in too_long[:5])
        raise ValueError(
            f"{genome}: {len(too_long)} sequence(s) longer than {MAX_SEQUENCE_LENGTH} bp: "
            f"{shown}. BAM cannot store these positions and minimap2 drops such "
            f"sequences; split them into pieces below 2^31 bp first.")


def build_hisat2_index(
    genome: Path,
    outdir: Path,
    threads: int = 4,
    prefix: str = "hisatidx",
) -> Path:
    """Build a HISAT2 index. Returns the index *prefix* path."""
    if not genome.is_file():
        raise FileNotFoundError(f"genome FASTA not found: {genome}")
    if shutil.which("hisat2-build") is None:
        raise RuntimeError("hisat2-build not found on PATH")

    outdir.mkdir(parents=True, exist_ok=True)
    idx_prefix = outdir / prefix

    cmd = [
        "hisat2-build",
        "-p", str(threads),
        str(genome),
        str(idx_prefix),
    ]
    log.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    log.info("HISAT2 index written with prefix %s", idx_prefix)
    return idx_prefix


def build_minimap2_index(
    genome: Path,
    outdir: Path,
    threads: int = 4,
    prefix: str = "mm2idx",
    part_bases: int | None = None,
) -> Path:
    """Build a minimap2 splice index. Returns the ``.mmi`` *file* path.

    Unlike HISAT2, minimap2 produces a single index file rather than a 6-file
    set, so the returned path is the ``.mmi`` itself rather than a stem.
    ``part_bases`` sets minimap2's ``-I`` (genome bases per index part;
    minimap2's default is 8 Gbp).
    """
    if not genome.is_file():
        raise FileNotFoundError(f"genome FASTA not found: {genome}")
    if shutil.which("minimap2") is None:
        raise RuntimeError("minimap2 not found on PATH")

    outdir.mkdir(parents=True, exist_ok=True)
    idx_path = outdir / f"{prefix}.mmi"

    cmd = [
        "minimap2",
        "-t", str(threads),
        "-x", "splice",
        *(["-I", str(int(part_bases))] if part_bases else []),
        "-d", str(idx_path),
        str(genome),
    ]
    log.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    log.info("minimap2 index written to %s", idx_path)
    parts = minimap2_index_parts(idx_path)
    if parts and len(parts) > 1:
        log.info("The index has %d parts (genome > ~8 Gbp); minimap2 loads one at a time "
                 "and VARUS aligns with --split-prefix.", len(parts))
    return idx_path

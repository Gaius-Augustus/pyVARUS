"""Download spots from SRA.

The online algorithm requires *spot-range* downloads (read N..X of run R), not
whole-run downloads. We therefore use:

* ``fastq-dump -N <n> -X <x> --fasta`` for batched range downloads (the proven
  legacy path);
* ``prefetch`` (``--prefetch`` mode) to fetch a run's ``.sra`` file once, after
  which ``fastq-dump -N/-X`` on the *local* file is fast (the remote
  spot-range path pays a multi-second resolver/HTTP latency per call and
  degrades for ranges deep inside large runs).
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

log = logging.getLogger(__name__)

# Number of retry attempts for transient SRA failures. The legacy code retries
# twice (3 attempts total) inside Downloader::getBatch.
DEFAULT_RETRIES = 3


@dataclass(frozen=True)
class BatchPaths:
    """Paths returned by ``download_batch``.

    For paired-end runs both ``r1`` and ``r2`` are populated. For single-end
    only ``r1`` is set; ``r2`` is ``None``.
    """
    r1: Path
    r2: Path | None
    batch_dir: Path

    def as_list(self) -> list[Path]:
        return [self.r1] if self.r2 is None else [self.r1, self.r2]


def batch_dir_for(outdir: Path, accession: str, n: int, x: int) -> Path:
    """Per-batch directory layout, identical to legacy ``Aligner::batchDir``::

        <outdir>/batches/<acc>/N<n>X<x>/
    """
    return outdir / "batches" / accession / f"N{n}X{x}"


def _require(tool: str) -> str:
    path = shutil.which(tool)
    if path is None:
        raise RuntimeError(f"{tool} not found on PATH")
    return path


def download_batch(
    accession: str,
    n: int,
    x: int,
    paired: bool,
    outdir: Path,
    *,
    retries: int = DEFAULT_RETRIES,
    fastq_dump: str = "fastq-dump",
    sra_path: Path | None = None,
) -> BatchPaths:
    """Download spots [n, x] of an SRA run as FASTA.

    Replicates ``legacy/Implementation/src/Downloader.cpp``'s
    ``shellCommand`` with ``fasta 120`` and ``--split-files`` for paired runs.

    ``sra_path`` (optional) points at a locally prefetched ``.sra`` file; it
    replaces the accession argument so ``fastq-dump`` reads from disk instead
    of the network. Output files are still named after the accession stem.

    Returns paths to the resulting FASTA file(s) inside the batch directory.
    """
    if shutil.which(fastq_dump) is None:
        raise RuntimeError(f"{fastq_dump} not found on PATH")

    bdir = batch_dir_for(outdir, accession, n, x)
    bdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        fastq_dump,
        "-N", str(n),
        "-X", str(x),
        "-O", str(bdir),
        "--fasta", "120",
    ]
    if paired:
        cmd.append("--split-files")
    cmd.append(str(sra_path) if sra_path is not None else accession)

    last_err: subprocess.CalledProcessError | None = None
    for attempt in range(1, retries + 1):
        log.info("fastq-dump %s N=%d X=%d (attempt %d/%d)%s",
                 accession, n, x, attempt, retries,
                 " [local .sra]" if sra_path is not None else "")
        try:
            subprocess.run(cmd, check=True)
            break
        except subprocess.CalledProcessError as e:
            last_err = e
            log.warning("fastq-dump failed (rc=%d) for %s N=%d X=%d",
                        e.returncode, accession, n, x)
    else:
        raise RuntimeError(
            f"fastq-dump failed after {retries} attempts for "
            f"{accession} N={n} X={x}: {last_err}"
        )

    if paired:
        r1 = bdir / f"{accession}_1.fasta"
        r2 = bdir / f"{accession}_2.fasta"
        # SRR097898 has a 3-file split (technical barcode in the middle); legacy
        # code handles that by using files[0] and files[-1]. We mirror that.
        if not r2.is_file():
            fastas = sorted(bdir.glob("*.fasta"))
            if len(fastas) >= 2:
                r1 = fastas[0]
                r2 = fastas[-1]
        if not r1.is_file() or not r2.is_file():
            raise RuntimeError(
                f"paired FASTA files missing in {bdir} for {accession}"
            )
        return BatchPaths(r1=r1, r2=r2, batch_dir=bdir)

    r1 = bdir / f"{accession}.fasta"
    if not r1.is_file():
        # A local .sra/.sralite input can yield a differently named file.
        fastas = sorted(bdir.glob("*.fasta"))
        if len(fastas) == 1:
            r1 = fastas[0]
    if not r1.is_file():
        raise RuntimeError(f"FASTA missing in {bdir} for {accession}")
    return BatchPaths(r1=r1, r2=None, batch_dir=bdir)


def find_prefetched(accession: str, sra_dir: Path) -> Path | None:
    """Locate a prefetched run file under ``sra_dir`` (``.sra`` preferred).

    sra-tools 3.x writes ``<sra_dir>/<ACC>/<ACC>.sra`` (or ``.sralite`` when
    the lite format is delivered). Older versions wrote ``<sra_dir>/<ACC>.sra``.
    """
    candidates = [
        sra_dir / accession / f"{accession}.sra",
        sra_dir / f"{accession}.sra",
        sra_dir / accession / f"{accession}.sralite",
        sra_dir / f"{accession}.sralite",
    ]
    for c in candidates:
        if c.is_file() and c.stat().st_size > 0:
            return c
    hits = sorted(sra_dir.glob(f"{accession}*/{accession}*.sra*")) + \
        sorted(sra_dir.glob(f"{accession}*.sra*"))
    for h in hits:
        if h.is_file() and h.stat().st_size > 0:
            return h
    return None


def prefetch_run(
    accession: str,
    sra_dir: Path,
    *,
    max_size_gb: float = 30.0,
    retries: int = DEFAULT_RETRIES,
    prefetch: str = "prefetch",
) -> Path:
    """Download the whole ``.sra`` of a run with ``prefetch``; return its path.

    ``--max-size`` bounds the transfer (prefetch refuses larger runs, which
    surfaces as a RuntimeError the controller treats as "keep using remote
    range dumps"). Existing complete files are reused.
    """
    _require(prefetch)
    sra_dir.mkdir(parents=True, exist_ok=True)
    existing = find_prefetched(accession, sra_dir)
    if existing is not None:
        log.info("prefetch %s: reusing %s", accession, existing)
        return existing

    cmd = [
        prefetch,
        "-O", str(sra_dir),
        "--max-size", f"{max(1, int(max_size_gb))}G",
        accession,
    ]
    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        log.info("prefetch %s (attempt %d/%d, max %sG)",
                 accession, attempt, retries, max(1, int(max_size_gb)))
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            last_err = e
            log.warning("prefetch failed (rc=%d) for %s", e.returncode, accession)
            continue
        found = find_prefetched(accession, sra_dir)
        if found is not None:
            return found
        last_err = RuntimeError("prefetch exited 0 but no .sra file was found")
    raise RuntimeError(
        f"prefetch failed after {retries} attempts for {accession}: {last_err}"
    )

"""``varus replay``: rebuild VARUS.bam from its manifest.

Reads ``VARUS.manifest.tsv`` and ``VARUS.splicedb.log.gz`` (see
:mod:`varus.provenance`), downloads the same spot ranges from SRA again,
aligns each batch against the splice-site DB version it was aligned with in
the original run, with the aligner thread count of the original run (HISAT2
2.2 places some reads differently with a different ``-p``), and merges the
batch BAMs into ``VARUS.bam``.

No sampling happens: the batches are fixed by the manifest. The result
holds the same alignments as the original BAM as long as SRA still serves
the runs and the same aligner version is used; the order of records with
equal coordinates may differ.
"""

from __future__ import annotations

import logging
import shutil
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from varus import __version__
from varus.align import (
    align_batch_hisat2,
    align_batch_minimap2,
    count_minimap2_quality,
    parse_hisat2_log,
)
from varus.download import download_batch
from varus.merge import merge_bams
from varus.provenance import (
    SPLICE_LOG_NAME,
    SpliceDBSnapshots,
    file_md5,
    read_manifest,
    tool_version,
)

log = logging.getLogger(__name__)

EXIT_INCOMPLETE = 2      # some batches could not be downloaded or aligned
MERGE_EVERY = 100        # batch BAMs per intermediate merge

# uniq_pct is written with 4 decimals; anything beyond rounding is a real
# difference (other reads, other aligner version, other DB).
_UNIQ_TOL = 0.01


@dataclass
class ReplayConfig:
    manifest: Path
    genome: Path
    index: Path
    outdir: Path
    threads: int = 4
    parallel_downloads: int = 4
    splice_db_log: Optional[Path] = None
    skip_genome_check: bool = False
    keep_batches: bool = False


def _check_genome(cfg: ReplayConfig, header: dict) -> None:
    want = header.get("genome_md5", "")
    if not want or cfg.skip_genome_check:
        if not want:
            log.warning("Manifest has no genome checksum; cannot verify %s", cfg.genome)
        return
    log.info("Checking the genome's MD5 (%s)", cfg.genome)
    got = file_md5(cfg.genome)
    if got != want:
        raise SystemExit(
            f"{cfg.genome}: MD5 {got or '(unreadable)'} differs from the genome the BAM "
            f"was built on ({header.get('genome', '?')}, MD5 {want}). Use the same "
            "assembly file; if it only differs in formatting (line width, compression "
            "removed), pass --skip-genome-check."
        )


def _check_versions(header: dict) -> None:
    if header.get("varus_version") and header["varus_version"] != __version__:
        log.warning("Manifest written by varus %s, replaying with %s",
                    header["varus_version"], __version__)
    aligner = header.get("aligner", "hisat2")
    for key, exe in (("aligner_version", aligner), ("samtools_version", "samtools")):
        want = header.get(key)
        have = tool_version(exe)
        if want and have != want:
            log.warning("%s version differs: manifest '%s', here '%s'; alignments may "
                        "differ slightly", exe, want, have)


def replay(cfg: ReplayConfig) -> int:
    """Rebuild VARUS.bam; returns 0, or EXIT_INCOMPLETE if batches are missing."""
    header, rows = read_manifest(cfg.manifest)
    if not rows:
        raise SystemExit(f"{cfg.manifest}: no batches; the original run wrote no BAM")
    longreads = header.get("mode") == "longreads"
    log_path = cfg.splice_db_log or (cfg.manifest.parent / header.get(
        "splice_db_log", SPLICE_LOG_NAME))
    need_db = any(r["db_version"] > 0 for r in rows)
    if need_db and not log_path.is_file():
        raise SystemExit(
            f"{log_path} not found. Batches were aligned with a splice-site DB whose "
            f"versions are recorded only there; pass it with --splice-db-log."
        )
    snaps = SpliceDBSnapshots(log_path) if need_db else None

    _check_genome(cfg, header)
    _check_versions(header)
    batch_size = int(header.get("batch_size") or 0)
    min_mapq = int(header.get("min_mapq") or (1 if longreads else 60))

    outdir = cfg.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    work = outdir / "replay"
    work.mkdir(exist_ok=True)
    K = max(1, int(cfg.parallel_downloads))
    log.info("Replaying %d downloads (%s batches) from %s", len(rows),
             header.get("batches", "?"), cfg.manifest)

    def fetch(r):
        return download_batch(r["accession"], r["n"], r["x"], r["paired"], work)

    bams: List[Path] = []
    parts: List[Path] = []
    missing: List[dict] = []
    n_diff = 0
    with ThreadPoolExecutor(max_workers=K, thread_name_prefix="replay-dl") as ex:
        pending: deque = deque()
        it = iter(rows)
        for r in it:
            pending.append((r, ex.submit(fetch, r)))
            if len(pending) >= K:
                break
        while pending:
            r, fut = pending.popleft()
            nxt = next(it, None)
            if nxt is not None:
                pending.append((nxt, ex.submit(fetch, nxt)))
            try:
                paths = fut.result()
            except Exception as e:
                log.error("batch %d (%s %d-%d): download failed: %s",
                          r["batch"], r["accession"], r["n"], r["x"], e)
                missing.append({**r, "reason": "download failed"})
                continue
            db = None
            if snaps is not None and r["db_version"] > 0:
                db = snaps.write(r["db_version"], paths.batch_dir / (
                    "intronDB.junc.bed" if longreads else "intronDB.splice_sites"))
            try:
                if longreads:
                    res = align_batch_minimap2(
                        reads=paths.r1, index=cfg.index, batch_dir=paths.batch_dir,
                        threads=r["align_threads"], preset=r["preset"] or "pacbio",
                        junc_bed=db, sort_threads=_sort_threads(cfg.threads),
                    )
                    uniq = count_minimap2_quality(res.bam, min_mapq=min_mapq)["uniq_pct"]
                else:
                    spots = batch_size if r["n_batches"] == 1 and batch_size else (
                        r["x"] - r["n"] + 1)
                    res = align_batch_hisat2(
                        r1=paths.r1, r2=paths.r2, index_prefix=cfg.index,
                        batch_dir=paths.batch_dir, threads=r["align_threads"], intron_db=db,
                        sort_threads=_sort_threads(cfg.threads),
                    )
                    uniq = parse_hisat2_log(res.log, batch_size=spots)["uniq_pct"]
            except RuntimeError as e:
                log.error("batch %d (%s %d-%d): alignment failed: %s",
                          r["batch"], r["accession"], r["n"], r["x"], e)
                missing.append({**r, "reason": "alignment failed"})
                continue
            if abs(uniq - r["uniq_pct"]) > _UNIQ_TOL:
                n_diff += 1
                log.warning("batch %d (%s %d-%d): uniquely mapped %.4f%%, original %.4f%%",
                            r["batch"], r["accession"], r["n"], r["x"], uniq, r["uniq_pct"])
            if not cfg.keep_batches:
                for p in paths.as_list() + ([db] if db else []):
                    Path(p).unlink(missing_ok=True)
            bams.append(res.bam)
            log.info("Replayed %d/%d: %s %d-%d", r["batch"], len(rows),
                     r["accession"], r["n"], r["x"])
            if len(bams) >= MERGE_EVERY:
                parts.append(merge_bams(bams, work / f"part_{len(parts):04d}.bam",
                                        threads=cfg.threads, compression=1))
                _remove(bams, cfg.keep_batches)
                bams = []

    inputs = parts + bams
    complete = not missing
    name = "VARUS.bam" if complete else "VARUS.incomplete.bam"
    if inputs:
        out = merge_bams(inputs, outdir / name, threads=cfg.threads)
        log.info("Replayed BAM: %s", out)
        _remove(inputs, cfg.keep_batches)
    if not cfg.keep_batches:
        shutil.rmtree(work / "batches", ignore_errors=True)
        _prune(work)
    if n_diff:
        log.warning("%d batch(es) differ in their uniquely-mapped share from the "
                    "original run (see warnings above)", n_diff)
    if not complete:
        miss = outdir / "replay_missing.tsv"
        with open(miss, "w", encoding="utf-8") as f:
            f.write("batch\taccession\tn\tx\treason\n")
            for r in missing:
                f.write(f"{r['batch']}\t{r['accession']}\t{r['n']}\t{r['x']}\t{r['reason']}\n")
        log.error("%d of %d downloads could not be replayed (listed in %s); the BAM "
                  "is incomplete and was written as %s", len(missing), len(rows), miss, name)
        return EXIT_INCOMPLETE
    return 0


def _sort_threads(threads: int) -> int:
    return max(1, min(4, threads - 1))


def _remove(paths: List[Path], keep: bool) -> None:
    if keep:
        return
    for p in paths:
        Path(p).unlink(missing_ok=True)


def _prune(d: Path) -> None:
    try:
        if d.is_dir() and not any(d.iterdir()):
            d.rmdir()
    except OSError:
        pass

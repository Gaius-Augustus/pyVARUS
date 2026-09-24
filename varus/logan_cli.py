"""argparse wiring for the ``varus logan`` subcommand.

Kept separate from :mod:`varus.cli` so the option table lives next to
:class:`varus.logan.LoganConfig`; ``cli.py`` only calls
:func:`add_logan_parser` and :func:`run_logan_cli`.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from varus.logan import TILE_WEIGHTS, LoganConfig, run_logan

log = logging.getLogger(__name__)


def add_logan_parser(sub: argparse._SubParsersAction) -> None:
    """Register the ``logan`` subcommand on an existing subparsers object."""
    p = sub.add_parser(
        "logan",
        help="Pre-screen and rank SRA runs using Logan contig assemblies.",
        description=(
            "Download Logan contigs for candidate runs from Runlist.tsv, "
            "spliced-align them to the genome, reject foreign/empty runs by "
            "tile breadth, greedily select complementary runs and write "
            "Runlist.logan.tsv plus intron hints."
        ),
    )
    p.add_argument("genome", type=Path, help="Genome FASTA file.")
    p.add_argument("--runlist", type=Path, required=True,
                   help="Runlist.tsv from 'varus runlist'.")
    p.add_argument("--outdir", type=Path, default=Path.cwd(),
                   help="Output directory; results go to <outdir>/logan/ and "
                        "<outdir>/Runlist.logan.tsv (default: cwd).")
    p.add_argument("--mmi", type=Path, default=None,
                   help="Existing minimap2 .mmi index (default: build one under "
                        "<outdir>/logan/genome/).")
    p.add_argument("--threads", type=int, default=4, help="minimap2/samtools threads.")
    p.add_argument("--download-workers", type=int, default=8,
                   help="Parallel HEAD/download connections to Logan S3.")
    p.add_argument("--max-candidates", type=int, default=500,
                   help="How many available runs to screen (0 = all).")
    p.add_argument("--chunk-runs", type=int, default=10,
                   help="Runs aligned per minimap2 invocation.")
    p.add_argument("--scan-workers", type=int, default=2,
                   help="Processes scanning chunk BAMs while the next chunk aligns "
                        "(0 = scan serially).")
    p.add_argument("--max-intron", type=int, default=20_000,
                   help="minimap2 -G maximum intron length.")
    p.add_argument("--min-contigs", type=int, default=100,
                   help="Runs with fewer contigs get status too_few_contigs.")
    p.add_argument("--min-tiles-frac", type=float, default=0.10,
                   help="Accept a run if it covers at least this fraction of the "
                        "tiles covered by the best run.")
    p.add_argument("--max-divergence", type=float, default=0.05,
                   help="Reject runs whose contigs' median divergence from the genome "
                        "(minimap2 de tag) exceeds this; catches other species whose "
                        "reads HISAT2 cannot map (0 = off).")
    p.add_argument("--ka-cap", type=float, default=50.0,
                   help="Cap on the per-contig k-mer abundance used as weight.")
    p.add_argument("--tile-weight", choices=list(TILE_WEIGHTS), default="ka_len",
                   help="Per-contig tile weight: unit, ka, or ka*len/150 (ka_len).")
    p.add_argument("--tile-size", type=int, default=5000, help="Tile size in bp.")
    p.add_argument("--select-top", type=int, default=50,
                   help="Maximum number of runs selected by the greedy ranking.")
    p.add_argument("--batch-size", type=int, default=50_000,
                   help="VARUS batch size used to scale the per-run pseudo-UMR mass.")
    p.add_argument("--prior-batches", type=float, default=1.0,
                   help="Pseudo-batches each run's contig evidence is worth.")
    p.add_argument("--keep-contigs", action="store_true",
                   help="Keep downloaded contig FASTAs under <outdir>/logan/contigs/.")
    p.add_argument("--logan-bam", action="store_true", dest="write_bam",
                   help="Keep chunk BAMs and write merged <outdir>/logan/LOGAN.bam.")
    p.add_argument("--seed", type=int, default=1, help="Random seed (backoff jitter).")
    p.add_argument("--timeout", type=float, default=60.0, help="HTTP timeout in seconds.")
    p.add_argument("--retries", type=int, default=4, help="HTTP retries per request.")
    p.add_argument("--longreads", action="store_true",
                   help="Record that the downstream 'varus run' stage uses long reads.")


def run_logan_cli(args: argparse.Namespace) -> int:
    """Build a :class:`LoganConfig` from parsed args and run the stage."""
    cfg = LoganConfig(
        genome=Path(args.genome),
        runlist=Path(args.runlist),
        outdir=Path(args.outdir),
        mmi=Path(args.mmi) if args.mmi else None,
        threads=args.threads,
        download_workers=args.download_workers,
        max_candidates=args.max_candidates,
        chunk_runs=args.chunk_runs,
        scan_workers=args.scan_workers,
        max_intron=args.max_intron,
        min_contigs=args.min_contigs,
        min_tiles_frac=args.min_tiles_frac,
        max_divergence=args.max_divergence,
        ka_cap=args.ka_cap,
        tile_weight=args.tile_weight,
        tile_size=args.tile_size,
        select_top=args.select_top,
        batch_size=args.batch_size,
        prior_batches=args.prior_batches,
        keep_contigs=args.keep_contigs,
        write_bam=args.write_bam,
        seed=args.seed,
        timeout=args.timeout,
        retries=args.retries,
        longreads=args.longreads,
    )
    if not cfg.genome.is_file():
        log.error("genome FASTA not found: %s", cfg.genome)
        return 2
    if not cfg.runlist.is_file():
        log.error("runlist not found: %s", cfg.runlist)
        return 2
    return run_logan(cfg)

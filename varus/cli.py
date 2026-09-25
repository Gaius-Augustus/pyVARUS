"""Command-line interface for VARUS.

Subcommands
-----------
runlist : query NCBI SRA for all RNA-seq runs of a species, write Runlist.tsv
index   : build a HISAT2 (default) or minimap2 (``--longreads``) index
logan   : pre-screen and rank runs by aligning their Logan contigs (optional)
run     : execute the online sampling loop (download + align + score)
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from varus import __version__


def _add_runlist(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "runlist",
        help="Fetch SRA RNA-seq run list for a species via NCBI Entrez.",
    )
    p.add_argument("species", help="Binomial species name, e.g. 'Drosophila melanogaster'")
    p.add_argument("--outdir", type=Path, default=Path.cwd(),
                   help="Output directory (default: cwd). Writes Runlist.tsv.")
    p.add_argument("--max-runs", type=int, default=0,
                   help="Limit to first N runs (0 = all available).")
    p.add_argument("--paired-only", action="store_true",
                   help="Keep only paired-end runs.")
    p.add_argument("--longreads", action="store_true",
                   help="Restrict the SRA query to long-read platforms "
                        "(PacBio SMRT, Oxford Nanopore). Without this flag the "
                        "runlist is dominated by Illumina short-read runs.")
    p.add_argument("--email", default=None,
                   help="Contact email for NCBI Entrez (recommended; "
                        "falls back to $NCBI_EMAIL).")
    p.add_argument("--api-key", default=None,
                   help="NCBI API key (optional; falls back to $NCBI_API_KEY).")


def _add_index(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "index",
        help="Build a HISAT2 (default) or minimap2 (--longreads) index.",
    )
    p.add_argument("genome", type=Path, help="Genome FASTA file.")
    p.add_argument("--outdir", type=Path, default=Path("genome"),
                   help="Output directory for the index (default: ./genome/).")
    p.add_argument("--threads", type=int, default=4,
                   help="Threads for the index builder (default: 4).")
    p.add_argument("--prefix", default=None,
                   help="Index file prefix (default: 'hisatidx' for short reads, "
                        "'mm2idx' for --longreads).")
    p.add_argument("--longreads", action="store_true",
                   help="Build a minimap2 splice index instead of HISAT2.")


def _add_run(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "run",
        help="Run the online sampling loop (download + align + score).",
    )
    p.add_argument("species",
                   help="Binomial species name (logged; the runs come from --runlist).")
    p.add_argument("genome", type=Path, help="Genome FASTA file.")
    p.add_argument("--runlist", type=Path, required=True, help="Path to Runlist.tsv.")
    p.add_argument("--index", type=Path, required=True,
                   help="HISAT2 index prefix (short reads, e.g. Sp/genome/hisatidx) "
                        "or minimap2 .mmi file (--longreads, e.g. Sp/genome/mm2idx.mmi).")
    p.add_argument("--outdir", type=Path, default=Path.cwd(),
                   help="Output directory.")
    p.add_argument("--batch-size", type=int, default=None,
                   help="Spots per batch (default: 50000 for short reads, "
                        "2000 for --longreads).")
    p.add_argument("--max-batches", type=int, default=1000)
    p.add_argument("--tile-size", type=int, default=5000)
    p.add_argument("--min-uniq-pct", type=float, default=5.0)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--keep-batches", action="store_true",
                   help="Keep per-batch FASTA/BAM files (default: delete after counting).")
    p.add_argument("--coverage-trace", type=int, default=0,
                   help="Write Coverage<N>.tsv every N batches (0 = never; final Coverage.csv "
                        "is always written).")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed (default: random).")
    p.add_argument("--bootstrap-all", action="store_true",
                   help="Download one batch from every run before starting the online loop "
                        "(equivalent to legacy --loadAllOnce).")
    p.add_argument("--profit-condition", action="store_true",
                   help="Stop early when expected profit ≤ 0. Off by default; matches the "
                        "legacy production setting (--profitCondition 0). The check is "
                        "always skipped on cold start (before any observations).")
    p.add_argument("--pipeline-downloads", action="store_true",
                   help="Keep one download in flight while the current batch is "
                        "aligned. Only meaningful with --parallel-downloads 1, "
                        "which otherwise reproduces the strictly serial v1 loop; "
                        "any --parallel-downloads K > 1 already pipelines.")
    p.add_argument("--longreads", action="store_true",
                   help="Align with minimap2 instead of HISAT2 (for PacBio Iso-Seq "
                        "or ONT direct-RNA). The platform is auto-detected per run "
                        "from the runlist (PACBIO_SMRT -> '-ax splice'; "
                        "OXFORD_NANOPORE -> '-ax splice -uf -k14'). Implies a "
                        "different splice-DB format and a smaller default --batch-size.")
    p.add_argument("--min-mapq", type=int, default=None,
                   help="MAPQ cutoff for the uniqueness gate (default: 60 for "
                        "short reads, 1 for --longreads).")
    p.add_argument("--advanced", nargs="*", default=[], metavar="KEY=VALUE",
                   help="Advanced overrides, e.g. lambda=3 pseudo-count=0.1 cost=0.001 "
                        "(defaults; v1 used lambda=10 pseudo-count=1).")

    # --- speed knobs (v2) ---
    g = p.add_argument_group("speed")
    g.add_argument("--parallel-downloads", type=int, default=6, metavar="K",
                   help="Keep K batch downloads in flight (default 6). K>1 implies "
                        "--pipeline-downloads; picks account for in-flight batches' "
                        "expected gains. K=1 reproduces the v1 pick sequence.")
    g.add_argument("--merge-batches", type=int, default=10, metavar="N",
                   help="Fetch up to N consecutive batches of a run in one download "
                        "(contiguous spot range) while greedy selection would pick "
                        "that run again anyway; only for runs whose first batch "
                        "passed the quality gate. Each fastq-dump call has a fixed "
                        "cost of 5-27 s. Default 10; 1 = off. Not used with "
                        "--parallel-downloads 1.")
    g.add_argument("--scan-workers", type=int, default=None, metavar="N",
                   help="Processes that scan a merged batch's BAM in parallel, split "
                        "by genome region (same counts as one pass). They get their "
                        "own share of --threads. Default: --threads/8, at most 4, "
                        "none below 16 threads; 0/1 = main thread.")
    g.add_argument("--no-align-ahead", action="store_true",
                   help="Do not align the next batch in the background while the "
                        "current one is scanned and scored (default: align ahead "
                        "whenever downloads are pipelined).")
    g.add_argument("--prefetch", action="store_true",
                   help="After a run has been picked --prefetch-after times, fetch its "
                        "whole .sra with `prefetch` in the background and range-dump "
                        "locally (removes the per-call remote latency).")
    g.add_argument("--prefetch-after", type=int, default=2)
    g.add_argument("--prefetch-max-gb", type=float, default=30.0,
                   help="Skip prefetch for runs estimated above this size (default 30).")
    g.add_argument("--prefetch-disk-gb", type=float, default=200.0,
                   help="Total disk budget for prefetched .sra files (default 200).")
    g.add_argument("--merge-every", type=int, default=100,
                   help="Merge batch BAMs in the background every N accepted batches "
                        "(default 100; 0 = single merge at the end).")
    g.add_argument("--no-hisat2-mm", action="store_true",
                   help="Do not pass --mm (memory-mapped index) to HISAT2.")
    g.add_argument("--keep-unaligned", action="store_true",
                   help="Keep unaligned reads in batch BAMs (default: hisat2 --no-unal).")
    g.add_argument("--splice-db-min-mult", type=int, default=1,
                   help="Only junctions with multiplicity >= N enter the aligner's "
                        "splice-site DB (default 1).")
    g.add_argument("--splice-db-rewrite-every", type=int, default=25,
                   help="Long-read BED12 DB refresh interval in batches (default 25).")

    # --- Logan pre-screen (v2) ---
    l = p.add_argument_group("logan")
    l.add_argument("--logan-dir", type=Path, default=None,
                   help="Output directory of `varus logan` (…/logan). Rejected runs are "
                        "dropped, the splice DB is seeded and accepted runs get an "
                        "estimator prior from their contig tile profile.")
    l.add_argument("--logan-top", type=int, default=0, metavar="K",
                   help="Keep only the K best-ranked Logan runs (0 = all accepted).")
    l.add_argument("--logan-only", action="store_true",
                   help="Also drop runs Logan could not process (absent/unsampled).")
    l.add_argument("--logan-prior-batches", type=float, default=1.0,
                   help="Weight of the Logan prior in batch equivalents (default 1; "
                        "0 = seed the splice DB only).")
    l.add_argument("--logan-prior-first-only", action="store_true",
                   help="Drop a run's Logan prior once it has a real batch, so a strong "
                        "prior (large --logan-prior-batches) only ranks unsampled runs.")
    l.add_argument("--no-logan-seed-db", action="store_true",
                   help="Do not seed intronDB from Logan introns.")
    l.add_argument("--logan-merge-introns", action="store_true",
                   help="Include Logan contig introns in the final introns.gff.")
    l.add_argument("--no-logan-bootstrap", action="store_true",
                   help="Do not give the first picks to the Logan-ranked runs in rank order.")
    l.add_argument("--logan-unprocessed-weight", type=float, default=-1.0,
                   help="Expected-read multiplier for runs Logan could not process "
                        "(default: the gate's acceptance rate; 1 = no discount).")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="varus",
        description="VARUS: online sampling of complementary RNA-seq reads from NCBI SRA.",
    )
    parser.add_argument("--version", action="version", version=f"varus {__version__}")
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    sub = parser.add_subparsers(dest="cmd", required=True)
    _add_runlist(sub)
    _add_index(sub)
    from varus.logan_cli import add_logan_parser
    add_logan_parser(sub)
    _add_run(sub)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=args.log_level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    if args.cmd == "runlist":
        from varus.runlist import fetch_runlist
        out = fetch_runlist(
            species=args.species,
            outdir=args.outdir,
            max_runs=args.max_runs,
            paired_only=args.paired_only,
            longreads=args.longreads,
            email=args.email,
            api_key=args.api_key,
        )
        print(out)
        return 0

    if args.cmd == "index":
        if args.longreads:
            from varus.index import build_minimap2_index
            out = build_minimap2_index(
                genome=args.genome,
                outdir=args.outdir,
                threads=args.threads,
                prefix=args.prefix or "mm2idx",
            )
        else:
            from varus.index import build_hisat2_index
            out = build_hisat2_index(
                genome=args.genome,
                outdir=args.outdir,
                threads=args.threads,
                prefix=args.prefix or "hisatidx",
            )
        print(out)
        return 0

    if args.cmd == "logan":
        from varus.logan_cli import run_logan_cli
        return run_logan_cli(args)

    if args.cmd == "run":
        from varus.controller import (
            Controller, VARUSConfig, apply_logan_prior, load_runs,
        )
        import random

        # Parse --advanced KEY=VALUE overrides
        advanced: dict[str, str] = {}
        for kv in (args.advanced or []):
            if "=" in kv:
                k, v = kv.split("=", 1)
                advanced[k.strip()] = v.strip()

        # Mode-dependent defaults: long-read SRA runs are smaller and the MAPQ
        # distribution from minimap2 is wider than HISAT2's.
        batch_size = args.batch_size
        if batch_size is None:
            batch_size = 2_000 if args.longreads else 50_000
        min_mapq = args.min_mapq
        if min_mapq is None:
            min_mapq = 1 if args.longreads else 60

        cfg = VARUSConfig(
            genome=args.genome,
            index_prefix=args.index,
            outdir=args.outdir,
            species=args.species,
            batch_size=batch_size,
            max_batches=args.max_batches,
            tile_size=args.tile_size,
            min_uniq_pct=args.min_uniq_pct,
            threads=args.threads,
            keep_batches=args.keep_batches,
            coverage_trace=args.coverage_trace,
            seed=args.seed,
            bootstrap_all=args.bootstrap_all,
            profit_condition=args.profit_condition,
            pipeline_downloads=args.pipeline_downloads,
            longreads=args.longreads,
            min_mapq=min_mapq,
            lambda_=float(advanced.get("lambda", 3.0)),
            pseudo_count=float(advanced.get("pseudo-count", 0.1)),
            cost=float(advanced.get("cost", 0.0)),
            parallel_downloads=max(1, args.parallel_downloads),
            align_ahead=not args.no_align_ahead,
            merge_batches=max(1, args.merge_batches),
            scan_workers=None if args.scan_workers is None else max(0, args.scan_workers),
            prefetch=args.prefetch,
            prefetch_after=args.prefetch_after,
            prefetch_max_gb=args.prefetch_max_gb,
            prefetch_disk_gb=args.prefetch_disk_gb,
            merge_every=args.merge_every,
            hisat2_mm=not args.no_hisat2_mm,
            keep_unaligned=args.keep_unaligned,
            splice_db_min_mult=args.splice_db_min_mult,
            splice_db_rewrite_every=args.splice_db_rewrite_every,
            logan_prior_batches=args.logan_prior_batches,
            logan_prior_first_only=args.logan_prior_first_only,
            logan_seed_db=not args.no_logan_seed_db,
            logan_merge_introns=args.logan_merge_introns,
            logan_bootstrap=not args.no_logan_bootstrap,
            logan_unprocessed_weight=args.logan_unprocessed_weight,
        )

        rng = random.Random(cfg.seed)
        runs = load_runs(args.runlist, cfg.batch_size, rng)
        if not runs:
            raise SystemExit("Runlist is empty or all runs are colorspace / filtered.")

        logan = None
        if args.logan_dir is not None:
            from varus.logan import load_logan
            logan = load_logan(args.logan_dir)
            runs = apply_logan_prior(
                runs, logan,
                batch_size=cfg.batch_size,
                prior_batches=cfg.logan_prior_batches,
                top=args.logan_top,
                only=args.logan_only,
            )
            if not runs:
                raise SystemExit(
                    "No runs left after applying the Logan pre-screen "
                    f"({args.logan_dir}); nothing to download."
                )
        ctrl = Controller(cfg, runs, logan=logan)
        return ctrl.run()

    raise SystemExit(f"unknown subcommand: {args.cmd}")


if __name__ == "__main__":
    sys.exit(main())

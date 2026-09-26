"""Command-line interface for VARUS.

Subcommands
-----------
runlist : query NCBI SRA for all RNA-seq runs of a species, write Runlist.tsv
index   : build a HISAT2 (default) or minimap2 (``--longreads``) index
logan   : pre-screen and rank runs by aligning their Logan contigs
run     : execute the online sampling loop (download + align + score);
          runs the Logan pre-screen first unless --no-logan or --logan-dir
replay  : rebuild VARUS.bam from VARUS.manifest.tsv + VARUS.splicedb.log.gz

``varus run`` and ``varus logan`` show only the options every user may need
in ``--help``; the expert options (sampling parameters, speed knobs, gates)
are listed by ``--help-all``.
"""

from __future__ import annotations

import argparse
import logging
import shlex
import shutil
import sys
from pathlib import Path
from typing import Optional

from varus import __version__


class HelpAllAction(argparse.Action):
    """``--help-all``: print the help including the expert options."""

    def __init__(self, option_strings, dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS, help=None):
        super().__init__(option_strings=option_strings, dest=dest,
                         default=default, nargs=0, help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        for a in parser._actions:
            h = getattr(a, "expert_help", None)
            if h is not None:
                a.help = h
        for g in parser._action_groups:
            d = getattr(g, "expert_description", None)
            if d is not None:
                g.description = d
        parser.print_help()
        parser.exit()


def expert_group(parser: argparse.ArgumentParser, title: str, description: str):
    """An argument group that ``--help`` hides and ``--help-all`` shows.

    The options work as usual; only their help text is suppressed until
    :class:`HelpAllAction` restores it.
    """
    g = parser.add_argument_group(title)
    g.expert_description = description  # type: ignore[attr-defined]
    orig_add = g.add_argument

    def add_argument(*args, **kwargs):
        a = orig_add(*args, **kwargs)
        a.expert_help = a.help or ""
        a.help = argparse.SUPPRESS
        return a

    g.add_argument = add_argument  # type: ignore[method-assign]
    return g


def add_help_all(parser: argparse.ArgumentParser, what: str) -> None:
    parser.add_argument("--help-all", action=HelpAllAction,
                        help=f"Show all options, including the expert options ({what}).")
    parser.epilog = (parser.epilog or "") + (
        f"Expert options ({what}) are listed by --help-all."
    )


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
        help="Run the Logan pre-screen and the online sampling loop.",
        description=(
            "Screen the candidate runs by their Logan contigs, then download, "
            "align and score read batches until --max-batches is reached. "
            "Writes VARUS.bam, introns.gff, Coverage.csv, RunStatistics.csv "
            "and BatchTimings.tsv to --outdir, plus VARUS.manifest.tsv and "
            "VARUS.splicedb.log.gz: archive these two (and the genome) to "
            "rebuild the BAM later with `varus replay`. Exits 3 when no batch "
            "passed the quality gate."
        ),
    )
    p.add_argument("species",
                   help="Binomial species name (logged; the runs come from --runlist).")
    p.add_argument("genome", type=Path, help="Genome FASTA file.")
    p.add_argument("--runlist", type=Path, required=True,
                   help="Runlist.tsv from 'varus runlist'.")
    p.add_argument("--index", type=Path, required=True,
                   help="HISAT2 index prefix (short reads, e.g. Sp/genome/hisatidx) "
                        "or minimap2 .mmi file (--longreads, e.g. Sp/genome/mm2idx.mmi).")
    p.add_argument("--outdir", type=Path, default=Path.cwd(),
                   help="Output directory (default: cwd).")
    p.add_argument("--threads", type=int, default=4,
                   help="CPU budget for everything that runs at once (default 4). "
                        "Set it to the cores the job owns; pyVARUS splits it "
                        "between the aligner, the scan and the downloads.")
    p.add_argument("--max-batches", type=int, default=1000,
                   help="Upper bound on downloaded batches (default 1000).")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed for a reproducible run order (default: random).")
    p.add_argument("--longreads", action="store_true",
                   help="Align with minimap2 instead of HISAT2 (PacBio Iso-Seq or "
                        "ONT direct-RNA runs). The preset is chosen per run from "
                        "the platform column of the runlist; --index is the .mmi "
                        "file; --batch-size defaults to 2000 and --min-mapq to 1.")
    p.add_argument("--no-logan", action="store_true",
                   help="Skip the Logan pre-screen. By default `varus run` runs it "
                        "before the first download (or reuses <outdir>/logan/ if "
                        "`varus logan` already wrote it); needs minimap2 on PATH.")
    add_help_all(p, "sampling parameters, speed knobs, Logan prior")

    g = expert_group(p, "expert: sampling", "Parameters of the online algorithm. "
                     "The defaults are those of the VARUS paper and the benchmarks.")
    g.add_argument("--batch-size", type=int, default=None,
                   help="Spots per batch (default: 50000 for short reads, "
                        "2000 for --longreads).")
    g.add_argument("--tile-size", type=int, default=5000,
                   help="Genome tile size in bp for the coverage score (default 5000).")
    g.add_argument("--min-uniq-pct", type=float, default=5.0,
                   help="Quality gate: reject a batch (and its run) when fewer than "
                        "this %% of its reads map uniquely (default 5).")
    g.add_argument("--min-mapq", type=int, default=None,
                   help="MAPQ cutoff for the uniqueness gate (default: 60 for "
                        "short reads, 1 for --longreads).")
    g.add_argument("--bootstrap-all", action="store_true",
                   help="Download one batch from every run before starting the online loop "
                        "(v1 --loadAllOnce).")
    g.add_argument("--profit-condition", action="store_true",
                   help="Stop early when the expected profit is <= 0 (v1 "
                        "--profitCondition 1). Off by default; never skipped a batch "
                        "in the benchmarks.")
    g.add_argument("--advanced", nargs="*", default=[], metavar="KEY=VALUE",
                   help="Estimator hyperparameters: lambda=3 pseudo-count=0.1 cost=0 "
                        "(defaults; v1 used lambda=10 pseudo-count=1).")
    g.add_argument("--keep-batches", action="store_true",
                   help="Keep per-batch FASTA/BAM files (default: delete after counting).")
    g.add_argument("--coverage-trace", type=int, default=0, metavar="N",
                   help="Write Coverage<N>.tsv every N batches (0 = never; the final "
                        "Coverage.csv is always written).")

    g = expert_group(p, "expert: speed", "Concurrency of downloads, alignment and "
                     "scanning. None of these changes what is sampled, except "
                     "--parallel-downloads 1, which reproduces the v1 pick sequence.")
    g.add_argument("--parallel-downloads", type=int, default=6, metavar="K",
                   help="Keep K batch downloads in flight (default 6); picks "
                        "account for in-flight batches' expected gains. K=1 is "
                        "the strictly serial v1 loop and pick sequence.")
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
    g.add_argument("--merge-every", type=int, default=100, metavar="N",
                   help="Merge batch BAMs in the background every N accepted batches "
                        "(default 100; 0 = single merge at the end).")
    g.add_argument("--splice-db-min-mult", type=int, default=1, metavar="N",
                   help="Only junctions seen >= N times enter the aligner's "
                        "splice-site DB (default 1).")

    g = expert_group(p, "expert: logan", "How the Logan pre-screen feeds the loop. "
                     "The pre-screen's own options (candidates, gates) belong to "
                     "`varus logan`; run it as a separate step to set them.")
    g.add_argument("--logan-dir", type=Path, default=None,
                   help="Use this `varus logan` output directory (.../logan) instead of "
                        "<outdir>/logan/. Rejected runs are dropped, the splice DB is "
                        "seeded and accepted runs get an estimator prior from their "
                        "contig tile profile.")
    g.add_argument("--logan-top", type=int, default=0, metavar="K",
                   help="Keep only the K best-ranked Logan runs (0 = all accepted).")
    g.add_argument("--logan-keep-unprocessed", action="store_true",
                   help="Keep the runs Logan could not process (newer than the last "
                        "Logan rebuild, or not among the screened candidates). By "
                        "default only accepted runs are sampled; the unprocessed ones "
                        "are kept anyway when the accepted runs hold fewer batches "
                        "than --max-batches (or when no run was accepted).")
    g.add_argument("--logan-prior-batches", type=float, default=1.0, metavar="B",
                   help="Weight of the Logan prior in batch equivalents (default 1; "
                        "0 = seed the splice DB only).")
    g.add_argument("--logan-merge-introns", action="store_true",
                   help="Include Logan contig introns in the final introns.gff.")


def _add_replay(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "replay",
        help="Rebuild VARUS.bam from its manifest (same reads, same alignment).",
        description=(
            "Download the batches listed in VARUS.manifest.tsv from SRA again and "
            "align each against the splice-site DB version it was aligned with "
            "(from VARUS.splicedb.log.gz, found next to the manifest). Writes "
            "VARUS.bam to --outdir; exits 2 and writes VARUS.incomplete.bam and "
            "replay_missing.tsv when batches could not be fetched."
        ),
    )
    p.add_argument("manifest", type=Path, help="VARUS.manifest.tsv of the original run.")
    p.add_argument("genome", type=Path,
                   help="Genome FASTA the original run used (checked by MD5).")
    p.add_argument("--index", type=Path, required=True,
                   help="Index of that genome from `varus index` (HISAT2 prefix, or the "
                        "minimap2 .mmi for a --longreads run).")
    p.add_argument("--outdir", type=Path, default=Path.cwd(),
                   help="Output directory (default: cwd).")
    p.add_argument("--threads", type=int, default=4,
                   help="Threads for sorting and merging (default 4). Each batch is "
                        "aligned with the thread count recorded in the manifest, "
                        "because HISAT2's output depends on it.")
    add_help_all(p, "downloads, file locations, checks")
    g = expert_group(p, "expert", "Rarely needed.")
    g.add_argument("--parallel-downloads", type=int, default=4, metavar="K",
                   help="Batch downloads in flight (default 4).")
    g.add_argument("--splice-db-log", type=Path, default=None,
                   help="VARUS.splicedb.log.gz if it is not next to the manifest.")
    g.add_argument("--skip-genome-check", action="store_true",
                   help="Do not compare the genome's MD5 with the manifest (for a "
                        "reformatted copy of the same assembly).")
    g.add_argument("--keep-batches", action="store_true",
                   help="Keep the per-batch FASTA/BAM files.")


def logan_prescreen(args: argparse.Namespace, cfg) -> Optional[Path]:
    """Run (or reuse) the Logan pre-screen for ``varus run``.

    Returns the Logan directory to load, or ``None`` when the loop should run
    without a prior: the pre-screen accepted no run (exit 3) or Logan's S3
    bucket was unreachable (exit 4). Any other failure is fatal.
    """
    log = logging.getLogger("varus.cli")
    ldir = Path(cfg.outdir) / "logan"
    if (ldir / "LoganRanking.tsv").is_file():
        log.info("Reusing the Logan pre-screen in %s", ldir)
        return ldir
    if shutil.which("minimap2") is None:
        raise SystemExit(
            "minimap2 not found on PATH; the Logan pre-screen needs it. "
            "Install minimap2 (conda install -c bioconda minimap2) or pass --no-logan."
        )
    from varus.logan import LoganConfig, run_logan
    lcfg = LoganConfig(
        genome=Path(cfg.genome),
        runlist=Path(args.runlist),
        outdir=Path(cfg.outdir),
        mmi=Path(cfg.index_prefix) if args.longreads else None,
        threads=cfg.threads,
        tile_size=cfg.tile_size,
        batch_size=cfg.batch_size,
        prior_batches=cfg.logan_prior_batches,
        seed=cfg.seed if cfg.seed is not None else 1,
    )
    log.info("Logan pre-screen: aligning candidate runs' contigs to %s", cfg.genome)
    rc = run_logan(lcfg)
    if rc == 0:
        return ldir
    if rc == 3:
        log.warning("Logan pre-screen accepted no run; sampling without a Logan prior "
                    "(rejected runs are still dropped)")
        return ldir if (ldir / "LoganRanking.tsv").is_file() else None
    if rc == 4:
        log.warning("Logan S3 unreachable; sampling without the pre-screen")
        return None
    raise SystemExit(f"varus logan failed with status {rc}")


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
    _add_replay(sub)
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

    if args.cmd in ("index", "run"):
        from varus.index import check_sequence_lengths
        check_sequence_lengths(args.genome)

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

    if args.cmd == "replay":
        from varus.replay import ReplayConfig, replay
        return replay(ReplayConfig(
            manifest=args.manifest,
            genome=args.genome,
            index=args.index,
            outdir=args.outdir,
            threads=max(1, args.threads),
            parallel_downloads=max(1, args.parallel_downloads),
            splice_db_log=args.splice_db_log,
            skip_genome_check=args.skip_genome_check,
            keep_batches=args.keep_batches,
        ))

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
            command=" ".join(shlex.quote(a) for a in ["varus"] + list(
                sys.argv[1:] if argv is None else argv)),
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
            longreads=args.longreads,
            min_mapq=min_mapq,
            lambda_=float(advanced.get("lambda", 3.0)),
            pseudo_count=float(advanced.get("pseudo-count", 0.1)),
            cost=float(advanced.get("cost", 0.0)),
            parallel_downloads=max(1, args.parallel_downloads),
            merge_batches=max(1, args.merge_batches),
            scan_workers=None if args.scan_workers is None else max(0, args.scan_workers),
            merge_every=args.merge_every,
            splice_db_min_mult=args.splice_db_min_mult,
            logan_prior_batches=args.logan_prior_batches,
            logan_merge_introns=args.logan_merge_introns,
        )

        rng = random.Random(cfg.seed)
        runs = load_runs(args.runlist, cfg.batch_size, rng)
        if not runs:
            raise SystemExit("Runlist is empty or all runs are colorspace / filtered.")

        logan = None
        logan_dir = args.logan_dir
        if logan_dir is None and not args.no_logan:
            logan_dir = logan_prescreen(args, cfg)
        if logan_dir is not None:
            from varus.logan import load_logan
            logan = load_logan(logan_dir)
            runs = apply_logan_prior(
                runs, logan,
                batch_size=cfg.batch_size,
                prior_batches=cfg.logan_prior_batches,
                top=args.logan_top,
                only=not args.logan_keep_unprocessed,
                max_batches=cfg.max_batches,
            )
            if not runs:
                raise SystemExit(
                    "No runs left after applying the Logan pre-screen "
                    f"({logan_dir}); nothing to download."
                )
            if not any(getattr(r, "logan_status", None) == "accepted" for r in runs):
                # Nothing accepted: keep the run filter, drop the prior.
                logan = None
        ctrl = Controller(cfg, runs, logan=logan)
        return ctrl.run()

    raise SystemExit(f"unknown subcommand: {args.cmd}")


if __name__ == "__main__":
    sys.exit(main())

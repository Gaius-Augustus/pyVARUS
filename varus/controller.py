"""VARUS online run-selection controller.

Implements the algorithm from Stanke et al. (2019) BMC Bioinformatics:
greedy maximization of the expected score gain

    S(c) = Σ_j ln(1 + c_j)

over 5 kb genome tiles j, where c_j is the UMR count in tile j. Each
iteration picks the SRA run whose next batch is expected to increase S the
most, downloads it, aligns it, and updates tile counts.

v2 additions (all off by default unless noted, see :class:`VARUSConfig`):

* ``parallel_downloads``: keep K batch downloads in flight. Picks are made
  by lazy greedy against ``total_obs`` *plus* the expected contribution of
  batches still in flight, so parallel picks are not blind repeats.
* ``prefetch``: fetch a run's ``.sra`` once and range-dump locally.
* Incremental splice-site DB: introns are stranded once and the aligner DB
  is rewritten only when new junctions appeared (on by default).
* Rolling merge of batch BAMs in the background (on by default).
* Logan prior: per-run pseudo-observations from ``varus logan`` inform the
  estimator before the first batch and seed the splice-site DB.
* ``run()`` returns an exit status; 3 = no batch passed the quality gate.

Key classes
-----------
VARUSConfig : All tunable parameters in one place.
RunState    : Per-run mutable state (observations, sigma, stats).
Controller  : Runs the online loop; owns global state (total_obs, introns).
"""

from __future__ import annotations

import logging
import math
import random
import shutil
import time
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

from varus.align import (
    align_batch_hisat2,
    align_batch_minimap2,
    count_minimap2_quality,
    parse_hisat2_log,
    preset_for_platform,
)
from varus.download import batch_dir_for, download_batch, find_prefetched, prefetch_run
from varus.estimator import AdvancedEstimator
from varus.introns import (
    IntronCounts,
    extract_introns_from_bam,
    read_introns_gff,
    write_introns_gff,
)
from varus.io import write_coverage, write_run_statistics
from varus.merge import merge_bams
from varus.runlist import RunRecord
from varus.strand import (
    StrandAssigner,
    assign_strand,
    write_hisat2_splice_sites,
    write_minimap2_junc_bed,
)
from varus.tiles import BAMStats, count_bam_stats

log = logging.getLogger(__name__)

Tile = Tuple[str, int]

# Exit status of ``varus run`` when every batch was rejected by the quality
# gate (or every download failed): no VARUS.bam can be written. Callers such
# as the BRAKER wrapper must treat this as "no usable RNA-seq", not a crash.
EXIT_NO_USABLE_DATA = 3

_TIMING_HEADER = (
    "batch\trun\tn\tx\tsuccess\tuniq_pct\tumrs\tt_download\tt_align\tt_scan\t"
    "t_db\tt_estimate\tt_batch\tinflight\tlocal_sra\n"
)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class VARUSConfig:
    """All controller parameters in one place."""
    genome: Path
    index_prefix: Path        # hisat2 index prefix (e.g. genome/hisatidx)
    outdir: Path

    batch_size: int = 50_000
    max_batches: int = 1_000
    tile_size: int = 5_000
    min_uniq_pct: float = 5.0  # quality gate; bad-quality if below

    threads: int = 4
    keep_batches: bool = False
    coverage_trace: int = 0    # write Coverage<N>.tsv every N batches; 0=never
    seed: Optional[int] = None
    bootstrap_all: bool = False

    # Estimator hyperparameters
    lambda_: float = 10.0
    pseudo_count: float = 1.0
    cost: float = 0.0          # per-read download cost (default 0 = ignore cost)
    # Stop early when expected profit ≤ 0. Off by default: matches the legacy
    # production pipeline (--profitCondition 0). When on, the check is also
    # skipped while no observations have been collected yet (cold start),
    # so the algorithm always gets at least one batch to bootstrap.
    profit_condition: bool = False

    # Pipeline downloads of round R+1 with alignments of round R. Downloads
    # are network-bound and single-threaded, alignments are CPU-bound — they
    # don't compete for the same resource. Equivalent to parallel_downloads=1
    # with one batch kept in flight. Default: off, so the algorithm matches
    # strict greedy ordering.
    pipeline_downloads: bool = False
    # Number of batch downloads kept in flight (>1 implies pipelining). Picks
    # for in-flight batches account for each other's expected tile gains.
    # NCBI tolerates a few concurrent fastq-dump streams; keep K <= 4.
    parallel_downloads: int = 1

    # Prefetch whole .sra files once a run has been picked `prefetch_after`
    # times, then range-dump locally. Bounded by per-run and total disk caps.
    prefetch: bool = False
    prefetch_after: int = 2
    prefetch_max_gb: float = 30.0
    prefetch_disk_gb: float = 200.0

    # Rolling merge: every `merge_every` accepted batches, merge their BAMs
    # into a part file in the background. 0 disables (single final merge).
    merge_every: int = 100

    # Aligner flags (see varus.align.align_batch_hisat2).
    hisat2_mm: bool = True
    keep_unaligned: bool = False

    # Splice-site DB maintenance. The DB is rewritten only when new stranded
    # junctions were added; long-read BED12 also carries multiplicity so it is
    # additionally refreshed every `splice_db_rewrite_every` batches.
    splice_db_rewrite_every: int = 25
    # Only junctions with multiplicity >= this enter the aligner DB.
    splice_db_min_mult: int = 1

    # Logan prior (populated by apply_logan_prior / the CLI).
    logan_prior_batches: float = 1.0
    logan_seed_db: bool = True
    logan_merge_introns: bool = False
    # First picks go to the Logan-ranked runs in rank order (one batch each)
    # before the estimator takes over. Without this the cold-start shared
    # prior (uniform over tiles) out-scores every informative profile.
    logan_bootstrap: bool = True
    # Expected-read multiplier for runs Logan could not process (absent,
    # unsampled). < 0 = use the gate's acceptance rate; 1.0 = no discount.
    logan_unprocessed_weight: float = -1.0

    # Long-read mode: align with minimap2 instead of HISAT2; feed back known
    # junctions via BED12 (--junc-bed) instead of HISAT2's tab format. The
    # minimap2 preset is chosen per run from the platform column of
    # Runlist.tsv (PACBIO_SMRT -> 'pacbio', OXFORD_NANOPORE -> 'ont').
    longreads: bool = False
    # Uniqueness MAPQ threshold for the quality gate. HISAT2 unique-mappers
    # all carry MAPQ=60, so this gate is effectively a no-op there; minimap2
    # emits a wider distribution, where MAPQ ≥ 1 excludes only ambiguous reads.
    min_mapq: int = 60


# ---------------------------------------------------------------------------
# Per-run state
# ---------------------------------------------------------------------------

@dataclass
class RunState:
    """Mutable state for a single SRA run."""
    record: RunRecord
    max_batches: int
    sigma: List[int]           # shuffled batch indices

    sigma_idx: int = 0
    times_downloaded: int = 0
    observations: Dict[Tile, int] = field(default_factory=dict)
    p: Dict[Tile, float] = field(default_factory=dict)
    expected_profit: float = 0.0
    avg_umr_pct: float = 0.0
    avg_spliced_pct: float = 0.0
    bad_quality: bool = False

    # v2: Logan pseudo-observations (tile -> pseudo-UMRs), array form of p
    # aligned to Controller._tiles, in-flight bookkeeping and prefetch state.
    prior_obs: Dict[Tile, float] = field(default_factory=dict)
    p_arr: Optional[np.ndarray] = field(default=None, repr=False)
    logan_status: str = ""          # accepted | unprocessed | "" (no Logan)
    logan_rank: int = 0
    logan_yield: Optional[float] = None   # 0..1, read-mass fraction on target
    n_inflight: int = 0
    sra_path: Optional[Path] = None
    prefetch_future: Optional[Future] = field(default=None, repr=False)
    prefetch_failed: bool = False
    last_pick: int = 0

    @classmethod
    def from_record(
        cls, record: RunRecord, batch_size: int, rng: random.Random
    ) -> "RunState":
        """Build a RunState from a RunRecord, initialising the sigma vector."""
        max_batches = max(1, math.ceil(record.total_spots / batch_size))
        sigma = list(range(max_batches))
        if len(sigma) > 1:
            # Shuffle all but the last element: mirrors shuffleExceptLast()
            # in legacy ChromosomeInitializer.cpp so the potentially-short
            # final batch is always downloaded last.
            front, tail = sigma[:-1], sigma[-1]
            rng.shuffle(front)
            sigma = front + [tail]
        return cls(record=record, max_batches=max_batches, sigma=sigma)

    @property
    def is_exhausted(self) -> bool:
        return self.sigma_idx >= len(self.sigma)

    def next_batch_range(self, batch_size: int) -> Tuple[int, int]:
        """Return (n, x) spot range for the next batch index in sigma."""
        k = self.sigma[self.sigma_idx]
        n = k * batch_size
        x = min((k + 1) * batch_size - 1, self.record.total_spots - 1)
        return n, x

    @property
    def estimated_sra_gb(self) -> float:
        """Rough size of the run's .sra file (2-bit bases + qualities)."""
        return self.record.total_bases * 1.1 / 4.0 / 1e9


# ---------------------------------------------------------------------------
# Per-batch task/result (returned from the parallel-safe download+align phase
# and consumed serially by _apply_batch_result)
# ---------------------------------------------------------------------------

@dataclass
class BatchTask:
    """One batch reserved for download.

    Holds the FASTA paths from a successful download. ``failed=True`` means
    the download itself failed (the run will be marked bad-quality when the
    task is consumed).
    """
    run: "RunState"
    n: int
    x: int
    paths: Optional["object"] = None      # download.BatchPaths
    failed: bool = False
    # v2 bookkeeping
    expected: Dict[Tile, float] = field(default_factory=dict)
    future: Optional[Future] = field(default=None, repr=False)
    seq: int = 0
    t_download: float = 0.0
    sra_path: Optional[Path] = None


@dataclass
class BatchResult:
    """Outcome of one download+align task. Mutated state lives in the run."""
    run: "RunState"
    success: bool
    bam_path: Optional[Path] = None
    bam_stats: Optional[BAMStats] = None
    introns: Optional[IntronCounts] = None
    uniq_pct: float = 0.0
    spliced_pct: float = 0.0
    t_align: float = 0.0
    t_scan: float = 0.0


# ---------------------------------------------------------------------------
# Controller
# ---------------------------------------------------------------------------

class Controller:
    """Online run-selection loop.

    Usage::

        cfg = VARUSConfig(genome=..., index_prefix=..., outdir=...)
        runs = [RunState.from_record(r, cfg.batch_size, rng) for r in records]
        ctrl = Controller(cfg, runs)
        status = ctrl.run()
    """

    def __init__(
        self,
        config: VARUSConfig,
        runs: List[RunState],
        logan: Optional["object"] = None,   # varus.logan.LoganPrior
    ) -> None:
        self.config = config
        self.runs = runs
        self.downloadable: List[RunState] = list(runs)
        self.rng = random.Random(config.seed)

        self.total_obs: Dict[Tile, int] = {}
        self.cumulative_introns: IntronCounts = IntronCounts({})
        # Stranded introns accumulated incrementally (feeds the aligner DB).
        self.stranded_introns: IntronCounts = IntronCounts({})
        self._batch_bams: List[Path] = []      # accepted batch BAMs not yet merged
        self._merged_parts: List[Path] = []
        self._merge_futures: List[Future] = []
        self._merge_ex: Optional[ThreadPoolExecutor] = None

        self.batch_count = 0
        self.total_score = 0.0
        self.total_profit = 0.0
        self.max_profit: float = 1.0  # initialised >0 so continuing() starts True

        # Running average of UMR% and spliced% across all downloaded batches.
        # Priors are 100% to overestimate initially (encourages exploration).
        self.avg_uniq: float = 100.0
        self.avg_spliced: float = 100.0

        self.estimator = AdvancedEstimator(
            lambda_=config.lambda_, pseudo_count=config.pseudo_count
        )

        # Tile index shared by the array paths of estimator and profit.
        self._tiles: List[Tile] = []
        self._tile_index: Dict[Tile, int] = {}
        self._x_arr: np.ndarray = np.zeros(0)
        self._logan_tiles: set = set()
        for r in runs:
            if r.prior_obs:
                self._logan_tiles.update(r.prior_obs.keys())

        # In-flight downloads and their simulated tile contributions.
        self._inflight: List[BatchTask] = []
        self._sim_extra: Dict[Tile, float] = {}
        self._seq = 0
        self._dl_ex: Optional[ThreadPoolExecutor] = None
        self._pf_ex: Optional[ThreadPoolExecutor] = None
        self._prefetched_bytes: Dict[str, int] = {}

        # Splice DB bookkeeping
        self._strander = StrandAssigner(config.genome)
        self._db_new_keys = 0
        self._db_batches_since_write = 0
        self._db_written = False

        self.config.outdir.mkdir(parents=True, exist_ok=True)
        self._splice_db_path = config.outdir / (
            "intronDB.junc.bed" if config.longreads else "intronDB.splice_sites"
        )
        self._timings_path = config.outdir / "BatchTimings.tsv"
        self._sra_dir = config.outdir / "sra"

        self._logan = logan
        self._logan_queue: List[RunState] = []
        self._unprocessed_weight = 1.0
        if logan is not None:
            self._seed_from_logan(logan)
            if config.logan_bootstrap:
                self._logan_queue = sorted(
                    (r for r in runs if r.logan_rank > 0), key=lambda r: r.logan_rank
                )
                if self._logan_queue:
                    log.info("Logan bootstrap: first %d picks follow the Logan ranking",
                             len(self._logan_queue))
            w = config.logan_unprocessed_weight
            if w < 0:
                rate = getattr(logan, "acceptance_rate", None)
                w = rate if rate is not None else 1.0
            self._unprocessed_weight = max(0.0, min(1.0, float(w)))
            if any(r.logan_status == "unprocessed" for r in runs):
                log.info("Logan: expected reads of unprocessed runs weighted by %.2f",
                         self._unprocessed_weight)

    # ------------------------------------------------------------------
    # Logan seed
    # ------------------------------------------------------------------

    def _seed_from_logan(self, logan) -> None:
        """Seed the splice DB and the strand cache from a ``varus logan`` run."""
        if not self.config.logan_seed_db:
            return
        introns_gff = getattr(logan, "introns", None)
        if introns_gff is not None and Path(introns_gff).is_file():
            seeded = read_introns_gff(Path(introns_gff))
            n_pre = self._strander.preload(seeded)
            self.stranded_introns = self.stranded_introns.merge(seeded)
            if self.config.logan_merge_introns:
                self.cumulative_introns = self.cumulative_introns.merge(seeded)
            log.info("Logan seed: %d stranded introns (%d cached)", len(seeded), n_pre)
            self._db_new_keys = len(seeded)
        src = getattr(logan, "junc_bed" if self.config.longreads else "splice_sites", None)
        if src is not None and Path(src).is_file():
            shutil.copyfile(Path(src), self._splice_db_path)
            self._db_written = True
            self._db_new_keys = 0
            log.info("Logan seed: splice DB copied to %s", self._splice_db_path)
        elif self.stranded_introns.counts:
            self._write_splice_db()

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def run(self) -> int:
        """Execute the online sampling loop. Returns an exit status."""
        if not self.downloadable:
            log.warning("No downloadable runs; nothing to do.")
            return self._finalize()

        if self.config.bootstrap_all:
            self._bootstrap()

        self._update_downloadable()
        self._estimate_p()
        self._calculate_profit()

        K = max(1, int(self.config.parallel_downloads))
        pipelined = self.config.pipeline_downloads or K > 1
        self._serial = not pipelined
        if not self._serial:
            self._dl_ex = ThreadPoolExecutor(
                max_workers=K, thread_name_prefix="varus-dl"
            )
        if self.config.prefetch:
            self._pf_ex = ThreadPoolExecutor(
                max_workers=1, thread_name_prefix="varus-prefetch"
            )
        if self.config.merge_every > 0:
            self._merge_ex = ThreadPoolExecutor(
                max_workers=1, thread_name_prefix="varus-merge"
            )
        self._max_inflight = 1 if self._serial else K
        if not self._timings_path.is_file():
            self._timings_path.write_text(_TIMING_HEADER, encoding="utf-8")

        try:
            self._refill()
            while self._inflight:
                t0 = time.monotonic()
                task = self._next_ready_task()
                self._inflight.remove(task)
                task.run.n_inflight -= 1
                self._release_expected(task)

                log.info(
                    "Batch %d/%d | downloadable=%d | inflight=%d | run: %s",
                    self.batch_count + 1,
                    self.config.max_batches,
                    len(self.downloadable),
                    len(self._inflight),
                    task.run.record.accession,
                )

                br = self._align_and_count(task)
                self.batch_count += 1
                self._apply_batch_result(br)

                t_db0 = time.monotonic()
                self._rebuild_intron_db()
                t_db = time.monotonic() - t_db0

                self.total_score = self._score()
                self.total_profit = (
                    self.total_score
                    - self.config.cost * self.config.batch_size * self.batch_count
                )

                t_est0 = time.monotonic()
                self._estimate_p()
                self._calculate_profit()
                t_est = time.monotonic() - t_est0
                self._export_stats(br.run)
                self._update_downloadable()
                self._write_timing(task, br, t_db, t_est, time.monotonic() - t0)

                if not self.downloadable:
                    log.info("No downloadable runs left.")
                    break
                if not self._continuing():
                    break
                self._refill()
        finally:
            self._drain_inflight()
            if self._dl_ex is not None:
                self._dl_ex.shutdown(wait=True)
            if self._pf_ex is not None:
                self._pf_ex.shutdown(wait=True)

        return self._finalize()

    # ------------------------------------------------------------------
    # Algorithm steps
    # ------------------------------------------------------------------

    def _bootstrap(self) -> None:
        """Download and process one batch from every run (--bootstrap-all)."""
        log.info("Bootstrap: processing one batch from each of %d runs", len(self.runs))
        for r in self.runs:
            if r.bad_quality or r.is_exhausted:
                continue
            n, x = r.next_batch_range(self.config.batch_size)
            r.sigma_idx += 1
            task = self._download_only(r, n, x)
            br = self._align_and_count(task)
            self.batch_count += 1
            self._apply_batch_result(br)
        self._rebuild_intron_db()
        self._estimate_p()

    def _estimate_p(self) -> None:
        """Re-estimate tile probability distributions for all runs."""
        tile_set = set(self.total_obs.keys())
        if self._logan_tiles:
            tile_set |= self._logan_tiles
        tiles = sorted(tile_set)
        if not tiles:
            return
        self._tiles = tiles
        self._tile_index = {t: i for i, t in enumerate(tiles)}
        self._x_arr = np.array(
            [self.total_obs.get(t, 0) for t in tiles], dtype=np.float64
        )
        arrays = self.estimator.estimate_arrays(
            tiles=tiles,
            obs_total=self.total_obs,
            run_obs=[r.observations for r in self.runs],
            times_downloaded=[r.times_downloaded for r in self.runs],
            prior_obs=[r.prior_obs or None for r in self.runs],
        )
        for run, arr in zip(self.runs, arrays):
            run.p_arr = arr
            # Dict form only for runs with real batches (cheap, keeps the
            # legacy attribute meaningful for tests and diagnostics).
            if run.times_downloaded > 0:
                run.p = dict(zip(tiles, arr.tolist()))

    def _calculate_profit(self) -> None:
        """Compute expectedProfit for every downloadable run; update avg stats."""
        n_stat = 4              # pseudocount (legacy numRuns=4)
        sum_umr = 100.0 * n_stat
        sum_spliced = 100.0 * n_stat

        no_reads_profit: Optional[float] = None

        for run in self.downloadable:
            shared = run.times_downloaded == 0 and not run.prior_obs
            if shared and no_reads_profit is not None:
                run.expected_profit = no_reads_profit
            else:
                run.expected_profit = self._profit(run)
                if shared:
                    no_reads_profit = run.expected_profit
                    log.debug(
                        "Prior profit (undownloaded runs): %.4f", no_reads_profit
                    )

            if run.times_downloaded > 0:
                sum_umr += run.avg_umr_pct
                sum_spliced += run.avg_spliced_pct
                n_stat += 1

        self.avg_uniq = sum_umr / n_stat
        self.avg_spliced = sum_spliced / n_stat

        if self.downloadable:
            self.max_profit = max(r.expected_profit for r in self.downloadable)

    def _effective_reads(self, run: RunState) -> float:
        """Expected useful reads of one batch (legacy: (umr% + spliced%) x B).

        For runs without batches the running averages are used; a Logan
        yield estimate scales them (mixed samples yield fewer reads), and
        runs Logan could not process are discounted by the gate's
        acceptance rate so unexplored runs do not crowd out known-good ones.
        """
        if run.times_downloaded > 0:
            umr_pct, spliced_pct = run.avg_umr_pct, run.avg_spliced_pct
            factor = 1.0
        else:
            umr_pct, spliced_pct = self.avg_uniq, self.avg_spliced
            if run.logan_yield is not None:
                factor = max(0.01, min(1.0, run.logan_yield))
            elif run.logan_status == "unprocessed":
                factor = self._unprocessed_weight
            else:
                factor = 1.0
        return factor * (umr_pct + spliced_pct) / 100.0 * self.config.batch_size

    def _profit(
        self, run: RunState, extra: Optional[Dict[Tile, float]] = None
    ) -> float:
        """Expected score gain from downloading one more batch of run r.

        Matches Controller::profit() in the legacy code. ``extra`` adds
        simulated observations (expected gains of in-flight batches) on top
        of ``total_obs`` so parallel picks account for each other.
        """
        effective = self._effective_reads(run)

        p_arr = run.p_arr
        if (
            p_arr is not None
            and self._tiles
            and p_arr.shape[0] == len(self._tiles)
        ):
            x = self._x_arr
            if extra:
                x = x.copy()
                for tile, v in extra.items():
                    i = self._tile_index.get(tile)
                    if i is not None:
                        x[i] += v
            pr = float(np.sum(np.log1p(x + p_arr * effective) - np.log1p(x)))
            return pr - self.config.cost * self.config.batch_size

        pr = 0.0
        for tile, prob in run.p.items():
            x = self.total_obs.get(tile, 0)
            if extra:
                x += extra.get(tile, 0.0)
            pr += math.log1p(x + prob * effective) - math.log1p(x)
        pr -= self.config.cost * self.config.batch_size
        return pr

    def _choose_next_run(self) -> Optional[RunState]:
        """Return the run with highest expectedProfit; break ties by avg_len.

        Ties are resolved by a weighted random draw proportional to avg_len
        (matching the legacy biasSelect() which prefers longer reads).
        """
        if not self.downloadable:
            return None

        best_profit = max(r.expected_profit for r in self.downloadable)
        self.max_profit = best_profit

        candidates = [
            r for r in self.downloadable if r.expected_profit == best_profit
        ]
        return self._tie_break(candidates)

    def _tie_break(self, candidates: List[RunState]) -> RunState:
        """Weighted random draw proportional to avg_len (legacy biasSelect)."""
        if len(candidates) == 1:
            return candidates[0]
        weights = [r.record.avg_len for r in candidates]
        total_w = sum(weights) or 1.0
        pick = self.rng.random() * total_w
        cumulative = 0.0
        for run, w in zip(candidates, weights):
            cumulative += w
            if pick <= cumulative:
                return run
        return candidates[-1]

    def _choose_next_run_simulated(self) -> Optional[RunState]:
        """Lazy-greedy pick accounting for in-flight batches.

        With nothing in flight this is exactly :meth:`_choose_next_run`.
        Otherwise candidates are visited in order of their (stale) expected
        profit; the profit of the current best is recomputed against
        ``total_obs + in-flight expectations`` and accepted once it is at
        least the next candidate's stale profit (Minoux's lazy greedy: the
        gain is submodular, so stale values are upper bounds).
        """
        if not self.downloadable:
            return None
        while self._logan_queue:
            run = self._logan_queue.pop(0)
            if run in self.downloadable and not run.is_exhausted and run.n_inflight == 0 \
                    and run.times_downloaded == 0:
                self.max_profit = max(self.max_profit, run.expected_profit)
                return run
        if not self._sim_extra:
            return self._choose_next_run()

        order = sorted(
            self.downloadable, key=lambda r: r.expected_profit, reverse=True
        )
        fresh: Dict[int, float] = {}
        # Runs with the shared prior have identical fresh profit; compute once.
        shared_fresh: Optional[float] = None
        i = 0
        while i < len(order):
            run = order[i]
            if run.times_downloaded == 0 and not run.prior_obs:
                if shared_fresh is None:
                    shared_fresh = self._profit(run, self._sim_extra)
                fresh[id(run)] = shared_fresh
            else:
                fresh[id(run)] = self._profit(run, self._sim_extra)
            next_stale = order[i + 1].expected_profit if i + 1 < len(order) else -math.inf
            if fresh[id(run)] >= next_stale - 1e-12:
                # Everything after `i` has stale <= fresh(run) -> run is optimal
                # among those; among already-recomputed ones pick the max,
                # breaking exact ties like the serial path does.
                return self._best_fresh(order[: i + 1], fresh)
            i += 1
        return self._best_fresh(order, fresh)

    def _best_fresh(self, cands: List[RunState], fresh: Dict[int, float]) -> RunState:
        best_val = max(fresh[id(r)] for r in cands)
        self.max_profit = best_val
        ties = [r for r in cands if abs(fresh[id(r)] - best_val) <= 1e-9]
        return self._tie_break(ties)

    def _continuing(self) -> bool:
        """True while the loop should keep running."""
        if self.config.max_batches > 0 and self.batch_count >= self.config.max_batches:
            log.info("Reached max_batches=%d; stopping.", self.config.max_batches)
            return False
        # Skip the profit check until we actually have observations. Without
        # this, the algorithm cannot bootstrap: every run starts with an empty
        # p distribution, profits are 0, and the loop would stop on iteration 1.
        if (
            self.config.profit_condition
            and self.total_obs
            and self.max_profit <= 0
        ):
            log.info("maxProfit=%.4f ≤ 0; stopping.", self.max_profit)
            return False
        return True

    def _score(self) -> float:
        """S(c) = Σ ln(1 + c_j) over all tiles with observations."""
        return sum(math.log1p(v) for v in self.total_obs.values())

    # ------------------------------------------------------------------
    # Download phase (network-bound; safe to overlap with alignment)
    # ------------------------------------------------------------------

    def _download_only(
        self, run: RunState, n: int, x: int, sra_path: Optional[Path] = None
    ) -> BatchTask:
        """Download one batch. No shared-state mutation."""
        t0 = time.monotonic()
        try:
            paths = download_batch(
                accession=run.record.accession,
                n=n, x=x,
                paired=run.record.paired,
                outdir=self.config.outdir,
                sra_path=sra_path,
            )
            task = BatchTask(run=run, n=n, x=x, paths=paths, failed=False)
        except RuntimeError as e:
            log.warning(
                "Download failed for %s [N=%d X=%d]: %s",
                run.record.accession, n, x, e,
            )
            task = BatchTask(run=run, n=n, x=x, paths=None, failed=True)
        task.t_download = time.monotonic() - t0
        task.sra_path = sra_path
        return task

    def _download_task(self, task: BatchTask) -> BatchTask:
        """Worker-thread body: fill ``task`` in place from the download."""
        done = self._download_only(task.run, task.n, task.x, task.sra_path)
        task.paths, task.failed = done.paths, done.failed
        task.t_download = done.t_download
        return task

    def _pick_and_download_single(self) -> Optional[BatchTask]:
        """Pick the best run, reserve its next sigma slot, and download one batch.

        Kept for the serial code path and for tests; the main loop uses
        :meth:`_refill`.
        """
        if not self.downloadable:
            return None
        if self.config.max_batches > 0 and self.batch_count >= self.config.max_batches:
            return None

        run = self._choose_next_run()
        if run is None or run.is_exhausted:
            return None

        n, x = run.next_batch_range(self.config.batch_size)
        run.sigma_idx += 1
        return self._download_only(run, n, x, self._local_sra(run))

    def _local_sra(self, run: RunState) -> Optional[Path]:
        """Resolve a finished prefetch into ``run.sra_path`` (main thread)."""
        fut = run.prefetch_future
        if fut is not None and fut.done():
            run.prefetch_future = None
            try:
                run.sra_path = fut.result()
                size = run.sra_path.stat().st_size if run.sra_path.is_file() else 0
                self._prefetched_bytes[run.record.accession] = size
                log.info("prefetch ready: %s (%.2f GB)",
                         run.sra_path, size / 1e9)
            except Exception as e:  # prefetch failed -> keep remote path
                run.prefetch_failed = True
                log.warning("prefetch failed for %s: %s (remote dumps continue)",
                            run.record.accession, e)
        if run.sra_path is not None and not run.sra_path.is_file():
            run.sra_path = None
        return run.sra_path

    def _maybe_prefetch(self, run: RunState) -> None:
        """Start a background prefetch when the trigger rule fires."""
        if not self.config.prefetch or self._pf_ex is None:
            return
        if run.sra_path is not None or run.prefetch_future is not None or run.prefetch_failed:
            return
        picks = run.times_downloaded + run.n_inflight
        if picks < self.config.prefetch_after:
            return
        est_gb = run.estimated_sra_gb
        if est_gb > self.config.prefetch_max_gb:
            log.info("prefetch skipped for %s: est. %.1f GB > cap %.1f GB",
                     run.record.accession, est_gb, self.config.prefetch_max_gb)
            run.prefetch_failed = True
            return
        used = sum(self._prefetched_bytes.values()) / 1e9
        if used + est_gb > self.config.prefetch_disk_gb:
            if not self._evict_prefetched(est_gb):
                log.info("prefetch deferred for %s: disk budget %.0f GB exhausted",
                         run.record.accession, self.config.prefetch_disk_gb)
                return
        existing = find_prefetched(run.record.accession, self._sra_dir) \
            if self._sra_dir.is_dir() else None
        if existing is not None:
            run.sra_path = existing
            self._prefetched_bytes[run.record.accession] = existing.stat().st_size
            return
        log.info("prefetch queued for %s (pick #%d, est. %.1f GB)",
                 run.record.accession, picks, est_gb)
        run.prefetch_future = self._pf_ex.submit(
            prefetch_run, run.record.accession, self._sra_dir,
            max_size_gb=self.config.prefetch_max_gb,
        )

    def _evict_prefetched(self, need_gb: float) -> bool:
        """Delete least-recently-picked prefetched runs until ``need_gb`` fits."""
        by_acc = {r.record.accession: r for r in self.runs}
        victims = sorted(
            (acc for acc in self._prefetched_bytes
             if acc in by_acc and by_acc[acc].n_inflight == 0),
            key=lambda acc: by_acc[acc].last_pick,
        )
        for acc in victims:
            used = sum(self._prefetched_bytes.values()) / 1e9
            if used + need_gb <= self.config.prefetch_disk_gb:
                return True
            self._drop_prefetched(by_acc[acc])
        used = sum(self._prefetched_bytes.values()) / 1e9
        return used + need_gb <= self.config.prefetch_disk_gb

    def _drop_prefetched(self, run: RunState) -> None:
        acc = run.record.accession
        self._prefetched_bytes.pop(acc, None)
        if run.sra_path is not None:
            d = run.sra_path.parent
            run.sra_path.unlink(missing_ok=True)
            if d != self._sra_dir and d.is_dir():
                shutil.rmtree(d, ignore_errors=True)
            run.sra_path = None
            log.info("prefetch evicted: %s", acc)

    def _refill(self) -> None:
        """Reserve picks and start downloads until K batches are in flight."""
        while len(self._inflight) < self._max_inflight:
            if not self.downloadable:
                return
            if (
                self.config.max_batches > 0
                and self.batch_count + len(self._inflight) >= self.config.max_batches
            ):
                return
            run = self._choose_next_run_simulated()
            if run is None or run.is_exhausted:
                self._update_downloadable()
                if run is None or not self.downloadable:
                    return
                continue

            n, x = run.next_batch_range(self.config.batch_size)
            run.sigma_idx += 1
            run.n_inflight += 1
            run.last_pick = self._seq
            self._seq += 1

            expected = self._expected_gain(run)
            for tile, v in expected.items():
                self._sim_extra[tile] = self._sim_extra.get(tile, 0.0) + v

            self._maybe_prefetch(run)
            sra_path = self._local_sra(run)

            if self._serial or self._dl_ex is None:
                task = self._download_only(run, n, x, sra_path)
                task.expected = expected
                task.seq = self._seq
            else:
                task = BatchTask(run=run, n=n, x=x, expected=expected,
                                 seq=self._seq, sra_path=sra_path)
                task.future = self._dl_ex.submit(self._download_task, task)
            self._inflight.append(task)
            self._update_downloadable()

    def _expected_gain(self, run: RunState) -> Dict[Tile, float]:
        """Expected per-tile UMRs of one batch of ``run`` (p × effective)."""
        eff = self._effective_reads(run)
        if run.p_arr is not None and self._tiles and run.p_arr.shape[0] == len(self._tiles):
            arr = run.p_arr * eff
            nz = np.nonzero(arr > 1e-9)[0]
            # Keep it sparse-ish: tiles carrying 99% of the mass or the top 5000.
            if nz.size > 5000:
                order = nz[np.argsort(arr[nz])[::-1][:5000]]
                nz = order
            return {self._tiles[i]: float(arr[i]) for i in nz}
        return {t: p * eff for t, p in run.p.items() if p * eff > 1e-9}

    def _release_expected(self, task: BatchTask) -> None:
        for tile, v in task.expected.items():
            cur = self._sim_extra.get(tile, 0.0) - v
            if cur <= 1e-9:
                self._sim_extra.pop(tile, None)
            else:
                self._sim_extra[tile] = cur

    def _next_ready_task(self) -> BatchTask:
        """Return the earliest-submitted task whose download has finished."""
        if self._serial or self._dl_ex is None:
            return self._inflight[0]
        futures = [t.future for t in self._inflight if t.future is not None]
        pending = [t for t in self._inflight if t.future is None]
        if pending:  # synchronously downloaded (should not happen in parallel mode)
            return pending[0]
        done, _ = wait(futures, return_when=FIRST_COMPLETED)
        ready = sorted(
            (t for t in self._inflight if t.future in done), key=lambda t: t.seq
        )
        task = ready[0]
        try:
            task.future.result()
        except Exception as e:  # defensive: _download_only catches its own errors
            log.warning("download worker raised for %s: %s", task.run.record.accession, e)
            task.failed = True
        return task

    def _drain_inflight(self) -> None:
        """Wait for (and discard) downloads still in flight at loop exit."""
        for t in self._inflight:
            if t.future is not None:
                try:
                    t.future.result()
                except Exception:
                    pass
            if t.paths is not None and not self.config.keep_batches:
                self._cleanup_batch_dir(t.paths, remove_bam=True)
        if self._inflight:
            log.info("Discarded %d in-flight download(s) at loop exit.", len(self._inflight))
        self._inflight = []
        self._sim_extra = {}

    # ------------------------------------------------------------------
    # Align + count phase (CPU-bound)
    # ------------------------------------------------------------------

    def _cleanup_batch_dir(self, paths, remove_bam: bool) -> None:
        """Delete FASTA and aligner logs (and optionally the BAM) of a batch.

        The aligner log is only needed for the quality gate, which has run by
        the time this is called. Empty directories are pruned so a 1000-batch
        run does not leave thousands of inodes behind.
        """
        for p in paths.as_list():
            Path(p).unlink(missing_ok=True)
        bdir = Path(paths.batch_dir)
        if bdir.is_dir():
            names = ["Log.final.out", "Log.minimap2.err"]
            if remove_bam:
                names.append("Aligned.out.bam")
            for name in names:
                (bdir / name).unlink(missing_ok=True)
        self._prune_dir(bdir)

    def _align_and_count(self, task: BatchTask) -> BatchResult:
        """Align one downloaded batch, count UMRs/introns, clean up FASTA."""
        run, n, x, paths = task.run, task.n, task.x, task.paths

        if task.failed or paths is None:
            if not self.config.keep_batches:
                bdir = batch_dir_for(self.config.outdir, run.record.accession, n, x)
                if bdir.is_dir():
                    shutil.rmtree(bdir, ignore_errors=True)
                    self._prune_dir(bdir.parent)
            return BatchResult(run=run, success=False)

        intron_db = (
            self._splice_db_path
            if self._splice_db_path.is_file()
            else None
        )
        threads = self.config.threads
        t0 = time.monotonic()
        try:
            if self.config.longreads:
                if paths.r2 is not None:
                    log.warning(
                        "Long-read run %s yielded paired FASTAs; "
                        "ignoring r2 and aligning r1 only.",
                        run.record.accession,
                    )
                preset = preset_for_platform(run.record.platform)
                result = align_batch_minimap2(
                    reads=paths.r1,
                    index=self.config.index_prefix,
                    batch_dir=paths.batch_dir,
                    threads=threads,
                    preset=preset,
                    junc_bed=intron_db,
                )
            else:
                result = align_batch_hisat2(
                    r1=paths.r1,
                    r2=paths.r2,
                    index_prefix=self.config.index_prefix,
                    batch_dir=paths.batch_dir,
                    threads=threads,
                    intron_db=intron_db,
                    mm=self.config.hisat2_mm,
                    keep_unaligned=self.config.keep_unaligned,
                )
        except RuntimeError as e:
            log.warning("Alignment failed for %s: %s", run.record.accession, e)
            if not self.config.keep_batches:
                self._cleanup_batch_dir(paths, remove_bam=True)
            return BatchResult(run=run, success=False, t_align=time.monotonic() - t0)
        t_align = time.monotonic() - t0

        t1 = time.monotonic()
        if self.config.longreads:
            stats = count_minimap2_quality(
                result.bam, min_mapq=self.config.min_mapq
            )
        else:
            stats = parse_hisat2_log(result.log, batch_size=self.config.batch_size)
        uniq_pct = stats["uniq_pct"]
        if uniq_pct < self.config.min_uniq_pct:
            log.warning(
                "Run %s batch [%d-%d]: uniq_pct=%.1f%% < min=%.1f%%; bad quality",
                run.record.accession, n, x, uniq_pct, self.config.min_uniq_pct,
            )
            if not self.config.keep_batches:
                self._cleanup_batch_dir(paths, remove_bam=True)
            return BatchResult(run=run, success=False, uniq_pct=uniq_pct,
                               t_align=t_align, t_scan=time.monotonic() - t1)

        bam_stats = count_bam_stats(result.bam, self.config.tile_size)
        spliced_pct = (
            100.0 * bam_stats.n_spliced / bam_stats.n_reads
            if bam_stats.n_reads > 0
            else 0.0
        )
        batch_introns = extract_introns_from_bam(result.bam)

        if not self.config.keep_batches:
            self._cleanup_batch_dir(paths, remove_bam=False)

        return BatchResult(
            run=run,
            success=True,
            bam_path=result.bam,
            bam_stats=bam_stats,
            introns=batch_introns,
            uniq_pct=uniq_pct,
            spliced_pct=spliced_pct,
            t_align=t_align,
            t_scan=time.monotonic() - t1,
        )

    def _apply_batch_result(self, br: BatchResult) -> None:
        """Serial state update from one BatchResult. Caller controls order."""
        run = br.run
        if not br.success:
            run.bad_quality = True
            if run.n_inflight == 0:
                self._drop_prefetched(run)
            return
        assert br.bam_stats is not None and br.introns is not None

        # Per-tile UMR accumulation
        for tile, count in br.bam_stats.umr_counts.items():
            run.observations[tile] = run.observations.get(tile, 0) + count
            self.total_obs[tile] = self.total_obs.get(tile, 0) + count

        run.times_downloaded += 1
        nd = run.times_downloaded
        run.avg_umr_pct += (br.uniq_pct - run.avg_umr_pct) / nd
        run.avg_spliced_pct += (br.spliced_pct - run.avg_spliced_pct) / nd

        log.info(
            "Run %s batch %d: uniq=%.1f%% spliced=%.1f%% UMRs=%d",
            run.record.accession, nd, br.uniq_pct, br.spliced_pct,
            sum(br.bam_stats.umr_counts.values()),
        )

        self.cumulative_introns = self.cumulative_introns.merge(br.introns)
        # Incremental stranding: only introns never seen before hit the genome.
        try:
            stranded, n_new = self._strander.assign_new(br.introns)
            self.stranded_introns = self.stranded_introns.merge(stranded)
            self._db_new_keys += n_new
        except Exception as e:
            log.warning("Incremental strand assignment failed: %s", e)
        self._db_batches_since_write += 1

        if br.bam_path is not None:
            self._batch_bams.append(br.bam_path)
            self._maybe_roll_merge()

    # ------------------------------------------------------------------
    # Splice-site DB
    # ------------------------------------------------------------------

    def _write_splice_db(self) -> int:
        introns = self.stranded_introns
        if self.config.splice_db_min_mult > 1:
            introns = IntronCounts({
                k: v for k, v in introns.counts.items()
                if v >= self.config.splice_db_min_mult
            })
        # The DB is only read by the aligner, which runs in the main thread
        # after this rewrite, so writing in place is race-free.
        if self.config.longreads:
            n = write_minimap2_junc_bed(introns, self._splice_db_path)
        else:
            n = write_hisat2_splice_sites(introns, self._splice_db_path)
        self._db_written = True
        self._db_new_keys = 0
        self._db_batches_since_write = 0
        return n

    def _rebuild_intron_db(self) -> None:
        """Refresh the aligner's splice-site DB from stranded introns.

        Incremental path: introns are stranded as batches arrive
        (:meth:`_apply_batch_result`); the DB file is rewritten only when
        new junctions appeared (HISAT2's tab format carries no multiplicity)
        or, for the BED12 long-read format, every ``splice_db_rewrite_every``
        batches. If ``stranded_introns`` is empty but ``cumulative_introns``
        is not (state set directly, e.g. on resume), fall back to a full
        rebuild via :func:`assign_strand`.

        Output format depends on the aligner: HISAT2 tab format for short reads,
        BED12 (``--junc-bed``) for minimap2 long reads.
        """
        if not self.cumulative_introns.counts and not self.stranded_introns.counts:
            return
        try:
            if not self.stranded_introns.counts and self.cumulative_introns.counts:
                stranded = assign_strand(self.cumulative_introns, self.config.genome)
                self.stranded_introns = stranded
                self._strander.preload(stranded)
                self._db_new_keys = len(stranded)
            need = self._db_new_keys > 0 or not self._db_written
            if self.config.longreads and self.config.splice_db_rewrite_every > 0:
                need = need or (
                    self._db_batches_since_write >= self.config.splice_db_rewrite_every
                )
            if need:
                n = self._write_splice_db()
                log.debug("Splice DB rewritten: %d junctions", n)
        except Exception as e:
            log.warning("Intron DB rebuild failed: %s", e)

    # ------------------------------------------------------------------
    # Rolling merge
    # ------------------------------------------------------------------

    def _maybe_roll_merge(self) -> None:
        if self._merge_ex is None or self.config.merge_every <= 0:
            return
        if len(self._batch_bams) < self.config.merge_every:
            return
        batch = list(self._batch_bams)
        self._batch_bams = []
        k = len(self._merge_futures)
        part = self.config.outdir / "merged" / f"part_{k:04d}.bam"
        self._merge_futures.append(
            self._merge_ex.submit(self._merge_part, batch, part)
        )
        log.info("Rolling merge #%d queued (%d BAMs)", k, len(batch))

    def _merge_part(self, bams: List[Path], part: Path) -> Path:
        merge_bams(
            bams, part,
            threads=max(1, self.config.threads // 2),
            compression=1,
        )
        if not self.config.keep_batches:
            for b in bams:
                b.unlink(missing_ok=True)
                self._prune_dir(b.parent)
        return part

    def _prune_dir(self, d: Path) -> None:
        """Remove ``d`` and its parents while empty, stopping at outdir."""
        stop = self.config.outdir.resolve()
        try:
            d = Path(d)
            while d.is_dir() and d.resolve() != stop and not any(d.iterdir()):
                d.rmdir()
                d = d.parent
        except OSError:
            pass

    # ------------------------------------------------------------------
    # Bookkeeping
    # ------------------------------------------------------------------

    def _update_downloadable(self) -> None:
        """Remove exhausted and bad-quality runs from the downloadable list."""
        before = len(self.downloadable)
        self.downloadable = [
            r for r in self.downloadable
            if not r.bad_quality and not r.is_exhausted
        ]
        removed = before - len(self.downloadable)
        if removed:
            log.info(
                "Removed %d run(s); %d remain downloadable.",
                removed, len(self.downloadable),
            )

    def _write_timing(
        self, task: BatchTask, br: BatchResult, t_db: float, t_est: float, t_batch: float
    ) -> None:
        umrs = sum(br.bam_stats.umr_counts.values()) if br.bam_stats else 0
        line = (
            f"{self.batch_count}\t{task.run.record.accession}\t{task.n}\t{task.x}\t"
            f"{int(br.success)}\t{br.uniq_pct:.2f}\t{umrs}\t{task.t_download:.2f}\t"
            f"{br.t_align:.2f}\t{br.t_scan:.2f}\t{t_db:.2f}\t{t_est:.2f}\t"
            f"{t_batch:.2f}\t{len(self._inflight)}\t{int(task.sra_path is not None)}\n"
        )
        log.info(
            "TIMING batch=%d run=%s dl=%.1fs align=%.1fs scan=%.1fs db=%.1fs est=%.1fs "
            "total=%.1fs inflight=%d",
            self.batch_count, task.run.record.accession, task.t_download,
            br.t_align, br.t_scan, t_db, t_est, t_batch, len(self._inflight),
        )
        try:
            with self._timings_path.open("a", encoding="utf-8") as f:
                f.write(line)
        except OSError as e:  # pragma: no cover
            log.debug("could not write timings: %s", e)

    # ------------------------------------------------------------------
    # Output helpers
    # ------------------------------------------------------------------

    def _export_stats(self, last_run: RunState) -> None:
        """Write per-batch coverage trace and run statistics."""
        cov_path = self.config.outdir / "Coverage.csv"
        stats_path = self.config.outdir / "RunStatistics.csv"

        # Coverage trace: every coverage_trace batches, write a snapshot
        if (
            self.config.coverage_trace > 0
            and self.batch_count % self.config.coverage_trace == 0
        ):
            trace = self.config.outdir / f"Coverage{self.batch_count}.tsv"
            write_coverage(self.total_obs, trace)

        # RunStatistics: every batch for the first 10, then every 10
        if self.batch_count < 10 or self.batch_count % 10 == 1:
            write_run_statistics(self.runs, stats_path)

    def _finalize(self) -> int:
        """Merge all batch BAMs, write final Coverage.csv and RunStatistics.csv.

        Returns the process exit status (0, or EXIT_NO_USABLE_DATA when no
        batch passed the quality gate).
        """
        # Wait for rolling merges.
        if self._merge_ex is not None:
            self._merge_ex.shutdown(wait=True)
            self._merge_ex = None
        for fut in self._merge_futures:
            try:
                self._merged_parts.append(fut.result())
            except Exception as e:
                log.error("Rolling merge failed: %s", e)
        self._merge_futures = []

        inputs = list(self._merged_parts) + list(self._batch_bams)
        log.info("Finalizing: merging %d BAM file(s) (%d parts + %d batches)",
                 len(inputs), len(self._merged_parts), len(self._batch_bams))

        write_coverage(self.total_obs, self.config.outdir / "Coverage.csv")
        write_run_statistics(self.runs, self.config.outdir / "RunStatistics.csv")

        status = 0
        if inputs:
            out_bam = self.config.outdir / "VARUS.bam"
            try:
                merge_bams(
                    inputs,
                    out_bam,
                    threads=self.config.threads,
                )
                log.info("Final BAM: %s", out_bam)
            except (RuntimeError, ValueError) as e:
                log.error("Final merge failed: %s", e)

            # Delete per-batch BAMs / parts unless the user asked to keep them
            if not self.config.keep_batches:
                for bam in inputs:
                    bam.unlink(missing_ok=True)
                    self._prune_dir(bam.parent)
        else:
            log.error(
                "No batch passed the quality gate (or every download failed); "
                "VARUS.bam not written. Exit status %d.", EXIT_NO_USABLE_DATA,
            )
            status = EXIT_NO_USABLE_DATA

        # Write the cumulative intron GFF alongside the final BAM
        if self.cumulative_introns.counts:
            gff_path = self.config.outdir / "introns.gff"
            write_introns_gff(self.cumulative_introns, gff_path)
            log.info("Cumulative introns: %s", gff_path)

        # Prefetched .sra files are scratch data.
        if not self.config.keep_batches and self._sra_dir.is_dir():
            shutil.rmtree(self._sra_dir, ignore_errors=True)
        self._strander.close()
        return status


# ---------------------------------------------------------------------------
# Runlist loader (called from CLI)
# ---------------------------------------------------------------------------

def load_runs(
    runlist_path: Path,
    batch_size: int,
    rng: random.Random,
    paired_only: bool = False,
) -> List[RunState]:
    """Parse a Runlist.tsv and return a list of RunState objects."""
    from varus.runlist import RunRecord

    records: List[RunRecord] = []
    with runlist_path.open(encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("@"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            acc = parts[0]
            try:
                spots = int(parts[1])
                bases = int(parts[2])
                avg_len = float(parts[3])
                paired = bool(int(parts[4]))
                colorspace = bool(int(parts[5]))
            except (ValueError, IndexError):
                log.warning("Skipping malformed runlist line: %s", line.rstrip())
                continue
            # 7th (platform) and 8th (bioproject) columns are optional for
            # backward compatibility with older Runlist.tsv files.
            platform = parts[6].strip() if len(parts) >= 7 else ""
            bioproject = parts[7].strip() if len(parts) >= 8 else ""
            if colorspace:
                continue
            if paired_only and not paired:
                continue
            records.append(
                RunRecord(
                    accession=acc,
                    total_spots=spots,
                    total_bases=bases,
                    avg_len=avg_len,
                    paired=paired,
                    colorspace=colorspace,
                    platform=platform,
                    bioproject=bioproject,
                )
            )

    log.info("Loaded %d runs from %s", len(records), runlist_path)
    return [RunState.from_record(r, batch_size, rng) for r in records]


# ---------------------------------------------------------------------------
# Logan prior (called from CLI)
# ---------------------------------------------------------------------------

def apply_logan_prior(
    runs: List[RunState],
    logan,                      # varus.logan.LoganPrior
    *,
    batch_size: int,
    prior_batches: float = 1.0,
    top: int = 0,
    only: bool = False,
) -> List[RunState]:
    """Filter and prime ``runs`` with the results of ``varus logan``.

    * Runs with status ``rejected`` are dropped.
    * ``top > 0`` keeps only the ``top`` best-ranked accepted runs (plus the
      never-processed ones unless ``only``).
    * Accepted runs get ``prior_obs``: their contig tile weights normalised
      to ``prior_batches × batch_size × 0.5`` pseudo-UMRs (0.5 ≈ typical
      UMR yield per read), i.e. worth ``prior_batches`` real batches.
    """
    status: Dict[str, str] = getattr(logan, "status", {}) or {}
    rank: Dict[str, int] = getattr(logan, "rank", {}) or {}
    tiles: Dict[str, Dict[Tile, float]] = getattr(logan, "tiles", {}) or {}
    yields: Dict[str, float] = getattr(logan, "yield_pct", {}) or {}
    scale = float(prior_batches) * batch_size * 0.5

    kept: List[RunState] = []
    n_rej = n_top = n_unproc = 0
    for r in runs:
        acc = r.record.accession
        st = status.get(acc, "unprocessed")
        if st == "rejected":
            n_rej += 1
            continue
        if st == "accepted":
            rk = rank.get(acc, 0)
            if top > 0 and (rk <= 0 or rk > top):
                n_top += 1
                continue
            r.logan_status = "accepted"
            r.logan_rank = rk
            y = yields.get(acc)
            if y is not None and y > 0:
                r.logan_yield = max(0.01, min(1.0, y / 100.0))
            t = tiles.get(acc)
            if t and scale > 0:
                total = float(sum(t.values()))
                if total > 0:
                    r.prior_obs = {k: v / total * scale for k, v in t.items()}
        else:
            n_unproc += 1
            if only:
                continue
            r.logan_status = "unprocessed"
        kept.append(r)
    log.info(
        "Logan prior: kept %d runs (%d rejected, %d outside top-%d, %d unprocessed%s)",
        len(kept), n_rej, n_top, top, n_unproc, " dropped" if only else " kept",
    )
    return kept

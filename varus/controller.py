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
import multiprocessing as mp
import os
import subprocess
import random
import shutil
import time
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from varus.align import (
    align_batch_hisat2,
    align_batch_minimap2,
    count_minimap2_quality,
    parse_hisat2_log,
    preset_for_platform,
    reserve_threads,
)
from varus.download import batch_dir_for, download_batch, find_prefetched, prefetch_run
from varus.estimator import AdvancedEstimator, sparse_counts
from varus.introns import (
    IntronCounts,
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
from varus.tiles import BAMStats, scan_batch_bam, scan_batch_bam_parallel

log = logging.getLogger(__name__)

Tile = Tuple[str, int]

# Exit status of ``varus run`` when every batch was rejected by the quality
# gate (or every download failed): no VARUS.bam can be written. Callers such
# as the BRAKER wrapper must treat this as "no usable RNA-seq", not a crash.
EXIT_NO_USABLE_DATA = 3

_TIMING_HEADER = (
    "batch\trun\tn\tx\tsuccess\tuniq_pct\tumrs\tt_download\tt_align\tt_scan\t"
    "t_db\tt_estimate\tt_batch\tinflight\tlocal_sra\tt_wait\tn_batches\n"
)


def auto_scan_workers(threads: int) -> int:
    """Default ``--scan-workers``: one scanner process per 8 threads, at most 4.

    48 threads give the benchmarked 4. Below 16 threads the scan stays in the
    main thread: there HISAT2 is the bottleneck and would lose the cores.
    """
    n = min(4, max(0, int(threads)) // 8)
    return n if n >= 2 else 0


def _sort_threads(align_threads: int) -> int:
    """``samtools sort -@`` for one batch: it sorts after the aligner has
    finished, and a 50 000-spot batch gains nothing beyond a few threads."""
    return max(1, min(4, align_threads - 1))


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class VARUSConfig:
    """All controller parameters in one place."""
    genome: Path
    index_prefix: Path        # hisat2 index prefix (e.g. genome/hisatidx)
    outdir: Path
    # Binomial species name from the CLI. Informational (logged at start);
    # the runs come from the runlist.
    species: str = ""

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
    # lambda 3 / a 0.1 (v1: 10 / 1): same or higher score on Drosophila and
    # both algae with and without Logan (docs/benchmark_logan.md).
    lambda_: float = 3.0
    pseudo_count: float = 0.1
    cost: float = 0.0          # per-read download cost (default 0 = ignore cost)
    # Stop early when expected profit ≤ 0. Off by default: matches the legacy
    # production pipeline (--profitCondition 0). When on, the check is also
    # skipped while no observations have been collected yet (cold start),
    # so the algorithm always gets at least one batch to bootstrap.
    profit_condition: bool = False

    # Number of batch downloads kept in flight (1 = strictly serial v1 loop). Picks
    # for in-flight batches account for each other's expected tile gains.
    # Default 6 (benchmark 2026-09: 1.46x over 3, 4x over 1); K=1 reproduces
    # the v1 pick sequence exactly.
    parallel_downloads: int = 6
    # Align the next downloaded batch in a background thread while the main
    # thread scans, applies and scores the current one. Only active with
    # pipelined downloads (the serial K=1 path must pick before downloading).
    # The aligner may then see a splice DB that lacks the current batch's
    # junctions (one batch stale); picks are unaffected.
    align_ahead: bool = True
    # Merge up to N consecutive picks of a run into one download of a
    # contiguous spot range (N × batch_size spots, one fastq-dump call, one
    # alignment). Extra batches are claimed only while greedy selection
    # would pick the same run again with the batches already claimed
    # counted as observed, and only for runs whose first batch passed the
    # quality gate. Each fastq-dump call has a fixed cost of 5-27 s, so runs
    # the loop exploits are fetched in a fraction of the time. 1 = off.
    # Default 10 (benchmark 2026-09-25, Drosophila: loop -37 % without Logan,
    # -63 % with Logan, same S). Never used with serial downloads, so K=1
    # still reproduces the v1 pick sequence.
    merge_batches: int = 10
    # Worker processes that scan a merged batch's BAM in parallel, split by
    # genome region (exact: reads spanning regions are merged centrally).
    # The scan is single-threaded Python and became the largest part of the
    # loop with merged batches. They get their own cores in the thread
    # budget. 0 or 1 = scan in the main thread; single batches always are.
    # None = auto_scan_workers(threads) (48 threads -> 4, < 16 -> 0).
    scan_workers: Optional[int] = None

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
    # Use the Logan prior only until a run's first batch; afterwards its own
    # observations replace it (the prior predicts the first batch).
    logan_prior_first_only: bool = False
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

    # Large runlists (mouse: 2 M runs): never-downloaded runs without a Logan
    # prior share one profit, so from this many of them on they are kept as
    # one pool (one representative for estimate/profit, expanded only in the
    # avg_len-weighted tie-break) instead of being scanned every batch.
    # 0 disables the pool.
    fresh_pool_min: int = 1000


# ---------------------------------------------------------------------------
# Per-run state
# ---------------------------------------------------------------------------

@dataclass
class RunState:
    """Mutable state for a single SRA run."""
    record: RunRecord
    max_batches: int
    # Shuffled batch indices; None until first needed when built lazily
    # (RunState.lazy), then drawn from Random(sigma_seed).
    sigma: Optional[List[int]]

    sigma_idx: int = 0
    times_downloaded: int = 0
    observations: Dict[Tile, int] = field(default_factory=dict)
    expected_profit: float = 0.0
    avg_umr_pct: float = 0.0
    avg_spliced_pct: float = 0.0
    bad_quality: bool = False

    # v2: Logan pseudo-observations (tile -> pseudo-UMRs), the run's tile
    # probability vector aligned to Controller._tiles, in-flight bookkeeping
    # and prefetch state.
    prior_obs: Dict[Tile, float] = field(default_factory=dict)
    p_arr: Optional[np.ndarray] = field(default=None, repr=False)
    # Sparse (indices, values) of observations + prior_obs against the
    # controller's tile index, rebuilt only when this run's counts change.
    obs_version: int = 0
    _sparse: Optional[tuple] = field(default=None, repr=False)
    logan_status: str = ""          # accepted | unprocessed | "" (no Logan)
    logan_rank: int = 0
    logan_yield: Optional[float] = None   # 0..1, read-mass fraction on target
    n_inflight: int = 0
    sra_path: Optional[Path] = None
    prefetch_future: Optional[Future] = field(default=None, repr=False)
    prefetch_failed: bool = False
    last_pick: int = 0
    sigma_seed: Optional[int] = None
    pos: int = -1                  # index in Controller.runs
    in_pool: bool = False          # member of the controller's fresh pool

    @staticmethod
    def _shuffled_sigma(max_batches: int, rng: random.Random) -> List[int]:
        sigma = list(range(max_batches))
        if len(sigma) > 1:
            # Shuffle all but the last element: mirrors shuffleExceptLast()
            # in legacy ChromosomeInitializer.cpp so the potentially-short
            # final batch is always downloaded last.
            front, tail = sigma[:-1], sigma[-1]
            rng.shuffle(front)
            sigma = front + [tail]
        return sigma

    @classmethod
    def from_record(
        cls, record: RunRecord, batch_size: int, rng: random.Random
    ) -> "RunState":
        """Build a RunState from a RunRecord, initialising the sigma vector."""
        max_batches = max(1, math.ceil(record.total_spots / batch_size))
        return cls(record=record, max_batches=max_batches,
                   sigma=cls._shuffled_sigma(max_batches, rng))

    @classmethod
    def lazy(cls, record: RunRecord, batch_size: int, seed: int) -> "RunState":
        """Like :meth:`from_record`, but the sigma vector is built on first use.

        A mouse runlist has 2 M runs and 860 M batch indices; building them
        all up front costs ~35 GB and minutes, while a run only ever picks
        from a few hundred of them.
        """
        max_batches = max(1, math.ceil(record.total_spots / batch_size))
        return cls(record=record, max_batches=max_batches, sigma=None,
                   sigma_seed=seed)

    def batch_order(self) -> List[int]:
        """The (shuffled) batch indices, built on first use."""
        if self.sigma is None:
            self.sigma = self._shuffled_sigma(
                self.max_batches, random.Random(self.sigma_seed))
        return self.sigma

    @property
    def is_exhausted(self) -> bool:
        if self.sigma is None:
            return self.sigma_idx >= self.max_batches
        return self.sigma_idx >= len(self.sigma)

    def next_batch_range(self, batch_size: int) -> Tuple[int, int]:
        """Return (n, x) spot range for the next batch index in sigma."""
        k = self.batch_order()[self.sigma_idx]
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
    # Alignment stage (possibly run ahead in the aligner thread).
    align_future: Optional[Future] = field(default=None, repr=False)
    aligned: Optional["_Aligned"] = None
    align_threads: int = 0     # set by the main thread from the thread budget
    n_batches: int = 1         # batches merged into this download (spot range)


@dataclass
class _Aligned:
    """Outcome of the aligner stage of one batch (no shared-state mutation)."""
    result: Optional[object] = None       # align.AlignmentResult
    t_align: float = 0.0
    error: str = ""


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
    n_batches: int = 1


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
        if config.scan_workers is None:
            config.scan_workers = auto_scan_workers(config.threads)
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
        self.max_profit: float = 1.0  # initialised >0 so continuing() starts True

        # Running average of UMR% and spliced% across all downloaded batches.
        # Priors are 100% to overestimate initially (encourages exploration).
        self.avg_uniq: float = 100.0
        self.avg_spliced: float = 100.0

        self.estimator = AdvancedEstimator(
            lambda_=config.lambda_, pseudo_count=config.pseudo_count
        )

        # Tile index shared by the array paths of estimator and profit. It
        # only grows (new tiles are appended in sorted order), so per-run
        # sparse counts stay valid between batches; ``_index_version`` is
        # bumped on the rare full rebuild, which invalidates them.
        self._tiles: List[Tile] = []
        self._tile_index: Dict[Tile, int] = {}
        self._index_version = 0
        self._x_arr: np.ndarray = np.zeros(0)
        self._log1p_x: np.ndarray = np.zeros(0)
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
        # Align-ahead: one batch aligning in the background while the main
        # thread processes the previous one. It keeps its in-flight
        # accounting (n_inflight, _sim_extra) until it is applied.
        self._align_ex: Optional[ThreadPoolExecutor] = None
        self._ahead: Optional[BatchTask] = None
        self._scan_ex = None   # ProcessPoolExecutor for merged-batch scans
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

        self._init_pool()

    # ------------------------------------------------------------------
    # Fresh pool: never-downloaded runs that share one profit
    # ------------------------------------------------------------------

    def _init_pool(self) -> None:
        """Move interchangeable fresh runs out of ``downloadable`` into a pool.

        A fresh run (no batch, no Logan prior, no yield estimate, not
        Logan-ranked) gets the shared prior from the estimator and hence the
        same expected profit as every other fresh run with the same Logan
        status. With 2 M such runs, scanning them in every estimate, profit
        and pick is what costs time, so they are represented by one member
        (``_pool_rep``) and expanded only in the avg_len-weighted tie-break,
        which draws exactly as the legacy tie-break over all of them would.
        """
        self._pool_n = 0
        self._pool_rep: Optional[RunState] = None
        self._pool_w: Optional[np.ndarray] = None
        self._pool_profit = -math.inf
        self._est_runs: List[RunState] = self.runs
        n_min = int(self.config.fresh_pool_min)
        if n_min <= 0 or self.config.bootstrap_all:
            return
        fresh = [r for r in self.downloadable if self._pool_eligible(r)]
        if len(fresh) < n_min:
            return
        counts: Dict[str, int] = {}
        for r in fresh:
            counts[r.logan_status] = counts.get(r.logan_status, 0) + 1
        status = max(counts, key=counts.get)
        for i, r in enumerate(self.runs):
            r.pos = i
        w = np.zeros(len(self.runs), dtype=np.float64)
        for r in fresh:
            if r.logan_status == status:
                r.in_pool = True
                w[r.pos] = r.record.avg_len
                self._pool_n += 1
        self._pool_w = w
        self._pool_rep = self.runs[self._first_pool_pos()]
        self.downloadable = [r for r in self.downloadable if not r.in_pool]
        self._est_runs = [r for r in self.runs if not r.in_pool]
        log.info("Fresh pool: %d interchangeable never-downloaded runs; %d runs tracked "
                 "individually", self._pool_n, len(self.downloadable))

    @staticmethod
    def _pool_eligible(r: RunState) -> bool:
        return (r.times_downloaded == 0 and not r.prior_obs and r.logan_yield is None
                and r.logan_rank <= 0 and not r.bad_quality and r.sigma_idx == 0
                and r.n_inflight == 0 and not r.is_exhausted)

    def _first_pool_pos(self) -> int:
        mask = self._pool_mask
        in_pool = np.flatnonzero(mask[self._pool_scan:])
        pos = self._pool_scan + int(in_pool[0])
        self._pool_scan = pos
        return pos

    @property
    def _pool_mask(self) -> np.ndarray:
        if getattr(self, "_pool_alive", None) is None:
            self._pool_alive = np.array([r.in_pool for r in self.runs], dtype=bool)
            self._pool_scan = 0
        return self._pool_alive

    def _pool_take(self, run: RunState) -> None:
        """``run`` was picked: track it individually from now on."""
        rep = self._pool_rep
        run.in_pool = False
        self._pool_alive[run.pos] = False
        self._pool_w[run.pos] = 0.0
        self._pool_n -= 1
        run.p_arr = rep.p_arr
        run.expected_profit = self._pool_profit
        # Keep downloadable in runs order (stable sorts and ties rely on it).
        i = len(self.downloadable)
        while i > 0 and self.downloadable[i - 1].pos > run.pos:
            i -= 1
        self.downloadable.insert(i, run)
        self._est_runs.append(run)
        if self._pool_n == 0:
            self._pool_rep = None
            self._pool_profit = -math.inf
        elif run is rep:
            new = self.runs[self._first_pool_pos()]
            new.p_arr = rep.p_arr
            new.expected_profit = rep.expected_profit
            self._pool_rep = new

    def _n_downloadable(self) -> int:
        return len(self.downloadable) + self._pool_n

    def _with_pool(self, runs: List[RunState]) -> List[RunState]:
        """``runs`` plus the pool representative (stands for the whole pool)."""
        return runs + [self._pool_rep] if self._pool_rep is not None else runs

    def _tie_break_pool(self, others: List[RunState]) -> RunState:
        """Tie-break over ``others`` plus every pool member, in runs order.

        Same draw as :meth:`_tie_break` on the full candidate list: one
        ``rng.random()`` scaled by the total avg_len, first candidate whose
        cumulative weight reaches it.
        """
        n = len(others) + self._pool_n
        if n == 1:
            return others[0] if others else self._pool_rep
        w = self._pool_w.copy()
        for r in others:
            w[r.pos] += r.record.avg_len
        cum = np.cumsum(w)
        total = float(cum[-1]) or 1.0
        pick = self.rng.random() * total
        i = int(np.searchsorted(cum, pick, side="left"))
        if i >= len(cum):
            last = max([r.pos for r in others] + [int(np.flatnonzero(self._pool_alive)[-1])])
            return self.runs[last]
        return self.runs[i]

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
        if self.config.species:
            log.info("Species: %s (%d runs in runlist)", self.config.species, len(self.runs))
        if not self._n_downloadable():
            log.warning("No downloadable runs; nothing to do.")
            return self._finalize()

        if self.config.bootstrap_all:
            self._bootstrap()

        self._update_downloadable()
        self._estimate_p()
        self._calculate_profit()

        K = max(1, int(self.config.parallel_downloads))
        self._serial = K == 1
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
        if not self._serial and self.config.align_ahead:
            self._align_ex = ThreadPoolExecutor(
                max_workers=1, thread_name_prefix="varus-align"
            )
        if self._parallel_scan_enabled():
            from concurrent.futures import ProcessPoolExecutor
            self._scan_ex = ProcessPoolExecutor(
                max_workers=int(self.config.scan_workers),
                mp_context=mp.get_context("spawn"),
            )
        self._max_inflight = 1 if self._serial else K
        self._log_thread_budget()
        if not self._timings_path.is_file():
            self._timings_path.write_text(_TIMING_HEADER, encoding="utf-8")

        try:
            self._refill()
            while self._inflight or self._ahead is not None:
                t0 = time.monotonic()
                if self._ahead is None:
                    self._start_align(self._take_ready_task())
                task, self._ahead = self._ahead, None
                self._await_align(task)
                t_wait = time.monotonic() - t0
                # Hand the aligner the next batch before processing this one.
                self._maybe_align_ahead()
                task.run.n_inflight -= 1
                self._release_expected(task)

                log.info(
                    "Batch %d/%d | downloadable=%d | inflight=%d | run: %s",
                    self.batch_count + 1,
                    self.config.max_batches,
                    self._n_downloadable(),
                    len(self._inflight),
                    task.run.record.accession,
                )

                br = self._count_aligned(task)
                br.n_batches = task.n_batches
                self.batch_count += task.n_batches
                self._apply_batch_result(br)

                t_db0 = time.monotonic()
                self._rebuild_intron_db()
                t_db = time.monotonic() - t_db0
                # A download may have finished meanwhile (the DB is fresh now).
                self._maybe_align_ahead()

                self.total_score = self._score()

                t_est0 = time.monotonic()
                self._estimate_p()
                self._calculate_profit()
                t_est = time.monotonic() - t_est0
                self._export_stats()
                self._update_downloadable()
                self._write_timing(task, br, t_db, t_est, time.monotonic() - t0, t_wait)

                if not self._n_downloadable():
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
            if self._align_ex is not None:
                self._align_ex.shutdown(wait=True)
                self._align_ex = None
            if self._scan_ex is not None:
                self._scan_ex.shutdown(wait=True)
                self._scan_ex = None

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

    def _sync_tile_index(self) -> bool:
        """Bring ``_tiles`` in line with the observed and Logan tiles.

        New tiles are appended in sorted order. Only when tiles disappeared
        (never in a real run; tests reset ``total_obs``) is the index rebuilt
        from scratch and ``_index_version`` bumped. Returns False when there
        are no tiles at all.
        """
        tile_set = set(self.total_obs.keys())
        if self._logan_tiles:
            tile_set |= self._logan_tiles
        if not tile_set:
            return False
        index = self._tile_index
        new = sorted(t for t in tile_set if t not in index)
        if len(self._tiles) + len(new) != len(tile_set):
            self._tiles = sorted(tile_set)
            self._tile_index = {t: i for i, t in enumerate(self._tiles)}
            self._index_version += 1
        elif new:
            base = len(self._tiles)
            self._tiles.extend(new)
            index.update((t, base + i) for i, t in enumerate(new))
        return True

    def _run_sparse(self, run: RunState) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """Sparse counts of ``run`` against the tile index (None = shared prior)."""
        prior = run.prior_obs
        if prior and run.times_downloaded > 0 and self.config.logan_prior_first_only:
            prior = {}
        if run.times_downloaded == 0 and not prior:
            return None
        key = (self._index_version, run.obs_version, id(run.prior_obs), len(run.prior_obs), bool(prior))
        cached = run._sparse
        if cached is not None and cached[0] == key:
            return cached[1]
        entry = sparse_counts(self._tile_index, run.observations, prior)
        run._sparse = (key, entry)
        return entry

    def _estimate_p(self) -> None:
        """Re-estimate tile probability distributions for all runs."""
        if not self._sync_tile_index():
            return
        tiles = self._tiles
        self._x_arr = np.array(
            [self.total_obs.get(t, 0) for t in tiles], dtype=np.float64
        )
        self._log1p_x = np.log1p(self._x_arr)
        runs = self._with_pool(self._est_runs)
        arrays = self.estimator.estimate_sparse(
            self._x_arr, [self._run_sparse(r) for r in runs]
        )
        for run, arr in zip(runs, arrays):
            run.p_arr = arr

    def _calculate_profit(self) -> None:
        """Compute expectedProfit for every downloadable run; update avg stats."""
        n_stat = 4              # pseudocount (legacy numRuns=4)
        sum_umr = 100.0 * n_stat
        sum_spliced = 100.0 * n_stat

        # Runs without counts share the estimator's prior array, so their
        # profit depends only on the expected reads (Logan yield/status).
        no_reads_profit: Dict[float, float] = {}

        def shared_profit(run: RunState) -> float:
            key = self._effective_reads(run)
            pr = no_reads_profit.get(key)
            if pr is None:
                pr = no_reads_profit[key] = self._profit(run)
                log.debug("Prior profit (undownloaded runs): %.4f", pr)
            return pr

        for run in self.downloadable:
            if run.times_downloaded == 0 and not run.prior_obs:
                run.expected_profit = shared_profit(run)
            else:
                run.expected_profit = self._profit(run)

            if run.times_downloaded > 0:
                sum_umr += run.avg_umr_pct
                sum_spliced += run.avg_spliced_pct
                n_stat += 1

        if self._pool_rep is not None:
            # Before the averages are updated, like the fresh runs above.
            self._pool_profit = shared_profit(self._pool_rep)
            self._pool_rep.expected_profit = self._pool_profit

        self.avg_uniq = sum_umr / n_stat
        self.avg_spliced = sum_spliced / n_stat

        runs = self._with_pool(self.downloadable)
        if runs:
            self.max_profit = max(r.expected_profit for r in runs)

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

    def _dense_extra(self, extra: Dict[Tile, float]) -> np.ndarray:
        """``extra`` as a vector over the tile index (unknown tiles dropped)."""
        arr = np.zeros(self._x_arr.shape[0], dtype=np.float64)
        index = self._tile_index
        for tile, v in extra.items():
            i = index.get(tile)
            if i is not None:
                arr[i] += v
        return arr

    def _profit(
        self,
        run: RunState,
        extra: Optional[Dict[Tile, float]] = None,
        x_extra: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    ) -> float:
        """Expected score gain from downloading one more batch of run r.

        Matches Controller::profit() in the legacy code. ``extra`` adds
        simulated observations (expected gains of in-flight batches) on top
        of ``total_obs`` so parallel picks account for each other;
        ``x_extra`` = ``(x, log1p(x))`` is the same thing already added to
        ``_x_arr`` (callers scoring many runs against one ``extra`` densify
        it once).
        """
        effective = self._effective_reads(run)

        p_arr = run.p_arr
        n = 0 if p_arr is None else p_arr.shape[0]
        if n and n <= self._x_arr.shape[0]:
            # p_arr may be shorter than the index when tiles were appended
            # since the last estimate; those tiles carry p = 0 for this run.
            if x_extra is not None:
                x, lx = x_extra[0][:n], x_extra[1][:n]
            elif extra:
                x = self._x_arr[:n] + self._dense_extra(extra)[:n]
                lx = np.log1p(x)
            else:
                x, lx = self._x_arr[:n], self._log1p_x[:n]
            pr = float(np.sum(np.log1p(x + p_arr * effective) - lx))
            return pr - self.config.cost * self.config.batch_size

        # No estimate yet (cold start before the first _estimate_p): the
        # expected gain is 0, only the cost applies.
        return -self.config.cost * self.config.batch_size

    def _choose_next_run(self) -> Optional[RunState]:
        """Return the run with highest expectedProfit; break ties by avg_len.

        Ties are resolved by a weighted random draw proportional to avg_len
        (matching the legacy biasSelect() which prefers longer reads).
        """
        if not self._n_downloadable():
            return None

        best_profit = max(r.expected_profit for r in self._with_pool(self.downloadable))
        self.max_profit = best_profit

        candidates = [
            r for r in self.downloadable if r.expected_profit == best_profit
        ]
        if self._pool_rep is not None and self._pool_profit == best_profit:
            return self._tie_break_pool(candidates)
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
        if not self._n_downloadable():
            return None
        while self._logan_queue:
            run = self._logan_queue.pop(0)
            if run in self.downloadable and not run.is_exhausted and run.n_inflight == 0 \
                    and run.times_downloaded == 0:
                self.max_profit = max(self.max_profit, run.expected_profit)
                return run
        if not self._sim_extra:
            return self._choose_next_run()

        # The pool representative stands for all pool members (same stale
        # and fresh profit); _best_fresh expands it. Sorting by position
        # among equal profits keeps the legacy stable order (rep is the
        # pool's first member).
        pool = self._pool_rep
        if pool is None:
            order = sorted(
                self.downloadable, key=lambda r: r.expected_profit, reverse=True
            )
        else:
            order = sorted(self._with_pool(self.downloadable),
                           key=lambda r: (-r.expected_profit, r.pos))
        fresh: Dict[int, float] = {}
        extra = self._sim_extra
        x_extra = None
        if self._x_arr.size:
            xe = self._x_arr + self._dense_extra(extra)
            x_extra = (xe, np.log1p(xe))
        # Runs with the shared prior have identical fresh profit per expected
        # read count; compute once.
        shared_fresh: Dict[float, float] = {}
        i = 0
        while i < len(order):
            run = order[i]
            if run.times_downloaded == 0 and not run.prior_obs:
                key = self._effective_reads(run)
                if key not in shared_fresh:
                    shared_fresh[key] = self._profit(run, extra, x_extra)
                fresh[id(run)] = shared_fresh[key]
            else:
                fresh[id(run)] = self._profit(run, extra, x_extra)
            next_stale = order[i + 1].expected_profit if i + 1 < len(order) else -math.inf
            at_pool = run is pool and self._pool_n > 1
            if at_pool:
                next_stale = run.expected_profit     # the next pool member
            if fresh[id(run)] >= next_stale - 1e-12:
                # Everything after `i` has stale <= fresh(run) -> run is optimal
                # among those; among already-recomputed ones pick the max,
                # breaking exact ties like the serial path does. Stopping at
                # the pool's first member means the other members were not
                # visited: rep stands for itself only.
                return self._best_fresh(order[: i + 1], fresh, expand=not at_pool)
            i += 1
        return self._best_fresh(order, fresh)

    def _best_fresh(
        self, cands: List[RunState], fresh: Dict[int, float], expand: bool = True
    ) -> RunState:
        best_val = max(fresh[id(r)] for r in cands)
        self.max_profit = best_val
        ties = [r for r in cands if abs(fresh[id(r)] - best_val) <= 1e-9]
        rep = self._pool_rep
        if expand and rep is not None and any(r is rep for r in ties):
            return self._tie_break_pool([r for r in ties if r is not rep])
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
            if not self._n_downloadable():
                return
            if (
                self.config.max_batches > 0
                and self.batch_count + self._pending_batches() >= self.config.max_batches
            ):
                return
            run = self._choose_next_run_simulated()
            if run is None or run.is_exhausted:
                self._update_downloadable()
                if run is None or not self._n_downloadable():
                    return
                continue
            if run.in_pool:
                self._pool_take(run)

            k0 = run.batch_order()[run.sigma_idx]
            n, x = run.next_batch_range(self.config.batch_size)
            run.sigma_idx += 1
            run.n_inflight += 1
            run.last_pick = self._seq
            self._seq += 1

            expected = self._expected_gain(run)
            m = self._merge_extra_batches(run, k0, expected)
            if m > 1:
                bs = self.config.batch_size
                x = min((k0 + m) * bs - 1, run.record.total_spots - 1)
                expected = {t: v * m for t, v in expected.items()}
            for tile, v in expected.items():
                self._sim_extra[tile] = self._sim_extra.get(tile, 0.0) + v

            self._maybe_prefetch(run)
            sra_path = self._local_sra(run)

            if self._serial or self._dl_ex is None:
                task = self._download_only(run, n, x, sra_path)
                task.expected = expected
                task.seq = self._seq
                task.n_batches = m
            else:
                task = BatchTask(run=run, n=n, x=x, expected=expected,
                                 seq=self._seq, sra_path=sra_path, n_batches=m)
                task.future = self._dl_ex.submit(self._download_task, task)
            self._inflight.append(task)
            self._update_downloadable()

    def _pending_batches(self) -> int:
        """Batches reserved but not yet applied (in flight + aligning ahead)."""
        n = sum(t.n_batches for t in self._inflight)
        if self._ahead is not None:
            n += self._ahead.n_batches
        return n

    def _merge_extra_batches(
        self, run: RunState, k0: int, expected: Dict[Tile, float]
    ) -> int:
        """How many batches (≥ 1) the pick of ``run`` at batch index ``k0`` covers.

        Claims batch indices k0+1, k0+2, ... (so the download stays one
        contiguous spot range) while they are still unused and greedy
        selection would pick ``run`` again with the claimed batches counted
        as observed (lazy greedy: the other runs' stale profits are upper
        bounds). Claimed indices are removed from the run's sigma.
        """
        m_max = int(self.config.merge_batches)
        if m_max <= 1 or run.times_downloaded < 1 or run.bad_quality:
            return 1
        if getattr(self, "_serial", True):
            return 1           # K=1 reproduces the v1 pick sequence
        budget = m_max
        if self.config.max_batches > 0:
            # this pick is already counted in neither batch_count nor pending
            budget = min(budget, self.config.max_batches - self.batch_count
                         - self._pending_batches())
        if budget <= 1 or not expected or not self._x_arr.size:
            return 1
        others = sorted(
            (r for r in self._with_pool(self.downloadable) if r is not run),
            key=lambda r: r.expected_profit, reverse=True,
        )
        xe = self._x_arr + self._dense_extra(self._sim_extra)
        step = self._dense_extra(expected)
        m = 1
        while m < budget:
            try:
                j = run.batch_order().index(k0 + m, run.sigma_idx)
            except ValueError:
                break          # next spot range already used: keep it contiguous
            xe = xe + step
            x_extra = (xe, np.log1p(xe))
            if not self._still_best(run, others, x_extra):
                break
            del run.sigma[j]
            m += 1
        if m > 1:
            log.debug("merged %d batches of %s from index %d", m, run.record.accession, k0)
        return m

    def _still_best(
        self, run: RunState, others: List[RunState], x_extra: Tuple[np.ndarray, np.ndarray]
    ) -> bool:
        """True if no other run beats ``run``'s fresh profit against ``x_extra``."""
        f_run = self._profit(run, None, x_extra)
        shared: Optional[float] = None
        for o in others:
            if o.expected_profit <= f_run + 1e-12:
                return True      # stale profits are upper bounds; the rest can't win
            if o.times_downloaded == 0 and not o.prior_obs:
                if shared is None:
                    shared = self._profit(o, None, x_extra)
                f_o = shared
            else:
                f_o = self._profit(o, None, x_extra)
            if f_o > f_run + 1e-12:
                return False
        return True

    def _expected_gain(self, run: RunState) -> Dict[Tile, float]:
        """Expected per-tile UMRs of one batch of ``run`` (p × effective)."""
        eff = self._effective_reads(run)
        if run.p_arr is not None and run.p_arr.size and run.p_arr.shape[0] <= len(self._tiles):
            arr = run.p_arr * eff
            nz = np.nonzero(arr > 1e-9)[0]
            # Keep it sparse-ish: tiles carrying 99% of the mass or the top 5000.
            if nz.size > 5000:
                order = nz[np.argsort(arr[nz])[::-1][:5000]]
                nz = order
            return {self._tiles[i]: float(arr[i]) for i in nz}
        return {}

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

    def _take_ready_task(self) -> BatchTask:
        """Block until a download is ready and remove it from the in-flight list."""
        task = self._next_ready_task()
        self._inflight.remove(task)
        return task

    def _maybe_align_ahead(self) -> None:
        """Start aligning the next finished download if the aligner is idle.

        Non-blocking: if no download has finished yet, the main loop picks
        the next batch up at the top of the next iteration.
        """
        if self._align_ex is None or self._ahead is not None or not self._inflight:
            return
        done = [t for t in self._inflight if t.future is not None and t.future.done()]
        if not done:
            return
        task = min(done, key=lambda t: t.seq)
        self._inflight.remove(task)
        try:
            task.future.result()
        except Exception as e:  # defensive: _download_only catches its own errors
            log.warning("download worker raised for %s: %s", task.run.record.accession, e)
            task.failed = True
        self._start_align(task)

    # ------------------------------------------------------------------
    # Thread budget: everything that runs at the same time shares --threads
    # ------------------------------------------------------------------

    @property
    def _merge_threads(self) -> int:
        """``samtools merge -@`` of a rolling merge (runs beside the loop)."""
        return max(1, min(4, self.config.threads // 8))

    def _merge_running(self) -> bool:
        return any(not f.done() for f in self._merge_futures)

    def _aligner_threads(self) -> int:
        """Threads for the next alignment, net of the work running beside it.

        Reserved: one core for the main thread while it scans and scores
        during an align-ahead alignment, one for the ``fastq-dump``
        processes when downloads are pipelined (they are latency-bound, so
        one core covers all K), and the rolling merge's threads while a
        merge runs. The final merge runs alone and gets all threads.
        """
        reserved = 0
        if self._align_ex is not None:
            # the main thread, or the scan workers it hands merged batches to
            reserved += max(1, self._n_scan_workers())
        if not getattr(self, "_serial", True):
            reserved += 1
        if self._merge_running():
            reserved += self._merge_threads
        return reserve_threads(self.config.threads, reserved)

    def _parallel_scan_enabled(self) -> bool:
        return (
            int(self.config.scan_workers) > 1
            and int(self.config.merge_batches) > 1
            and not self._serial
            and not self.config.longreads   # the exact split relies on HISAT2's NH tags
        )

    def _n_scan_workers(self) -> int:
        return int(self.config.scan_workers) if self._scan_ex is not None else 0

    def _log_thread_budget(self) -> None:
        a = self._aligner_threads()
        nsw = self._n_scan_workers()
        log.info(
            "Thread budget (--threads %d): aligner %d (samtools sort -@ %d)%s%s; "
            "rolling merge -@ %d, taken from the aligner while it runs",
            self.config.threads, a, _sort_threads(a),
            (f", scan workers {nsw} (merged batches; main thread otherwise)" if nsw
             else ", main thread 1 (align-ahead)") if self._align_ex is not None else "",
            f", downloads 1 ({self._max_inflight} fastq-dump)" if not self._serial else "",
            self._merge_threads,
        )

    def _start_align(self, task: BatchTask) -> None:
        """Align ``task`` in the aligner thread (or inline without one)."""
        task.align_threads = self._aligner_threads()
        if self._align_ex is not None:
            task.align_future = self._align_ex.submit(self._align_only, task)
        else:
            task.aligned = self._align_only(task)
        self._ahead = task

    def _await_align(self, task: BatchTask) -> _Aligned:
        if task.align_future is not None:
            try:
                task.aligned = task.align_future.result()
            except Exception as e:  # defensive: _align_only catches RuntimeError
                task.aligned = _Aligned(error=str(e))
            task.align_future = None
        if task.aligned is None:
            task.aligned = self._align_only(task)
        return task.aligned

    def _drain_inflight(self) -> None:
        """Wait for (and discard) downloads still in flight at loop exit."""
        pending = list(self._inflight)
        if self._ahead is not None:
            pending.append(self._ahead)
            if self._ahead.align_future is not None:
                try:
                    self._ahead.align_future.result()
                except Exception:
                    pass
            self._ahead = None
        for t in pending:
            if t.future is not None:
                try:
                    t.future.result()
                except Exception:
                    pass
            if t.paths is not None and not self.config.keep_batches:
                self._cleanup_batch_dir(t.paths, remove_bam=True)
        if pending:
            log.info("Discarded %d in-flight batch(es) at loop exit.", len(pending))
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
        task.aligned = self._align_only(task)
        return self._count_aligned(task)

    def _align_only(self, task: BatchTask) -> _Aligned:
        """Run the aligner on one downloaded batch.

        Safe to run in the aligner thread: it reads only the task, the config
        and the splice-DB file (replaced atomically by the main thread) and
        leaves all cleanup and state updates to :meth:`_count_aligned`.
        """
        run, paths = task.run, task.paths
        if task.failed or paths is None:
            return _Aligned(error="download failed")

        intron_db = (
            self._splice_db_path
            if self._splice_db_path.is_file()
            else None
        )
        threads = task.align_threads or self.config.threads
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
                    sort_threads=_sort_threads(threads),
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
                    sort_threads=_sort_threads(threads),
                )
        except RuntimeError as e:
            return _Aligned(t_align=time.monotonic() - t0, error=str(e) or "alignment failed")
        if self._scan_ex is not None and task.n_batches > 1:
            # index for the region-split scan (here, not in the main thread);
            # CSI because BAI cannot hold chromosomes over 512 Mbp (wheat 3B)
            try:
                subprocess.run(["samtools", "index", "-c", "-@", "2", str(result.bam)],
                               check=True, capture_output=True)
            except (OSError, subprocess.CalledProcessError) as e:
                log.debug("samtools index failed (%s); the scan indexes itself", e)
        return _Aligned(result=result, t_align=time.monotonic() - t0)

    def _count_aligned(self, task: BatchTask) -> BatchResult:
        """Quality gate, BAM scan and cleanup of an aligned batch (main thread)."""
        run, n, x, paths = task.run, task.n, task.x, task.paths
        aligned = task.aligned if task.aligned is not None else self._await_align(task)

        if task.failed or paths is None:
            if not self.config.keep_batches:
                bdir = batch_dir_for(self.config.outdir, run.record.accession, n, x)
                if bdir.is_dir():
                    shutil.rmtree(bdir, ignore_errors=True)
                    self._prune_dir(bdir.parent)
            return BatchResult(run=run, success=False)

        t_align = aligned.t_align
        if aligned.result is None:
            log.warning("Alignment failed for %s: %s", run.record.accession, aligned.error)
            if not self.config.keep_batches:
                self._cleanup_batch_dir(paths, remove_bam=True)
            return BatchResult(run=run, success=False, t_align=t_align)
        result = aligned.result

        t1 = time.monotonic()
        if self.config.longreads:
            stats = count_minimap2_quality(
                result.bam, min_mapq=self.config.min_mapq
            )
        else:
            spots = self.config.batch_size if task.n_batches == 1 else (x - n + 1)
            stats = parse_hisat2_log(result.log, batch_size=spots)
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

        if self._scan_ex is not None and task.n_batches > 1:
            try:
                bam_stats, batch_introns = scan_batch_bam_parallel(
                    result.bam, self.config.tile_size, self._scan_ex,
                    n_parts=4 * self._n_scan_workers(),
                )
            except Exception as e:  # e.g. a dead worker: fall back to one pass
                log.warning("parallel scan failed (%s); scanning in the main thread", e)
                bam_stats, batch_introns = scan_batch_bam(result.bam, self.config.tile_size)
            for ext in (".bai", ".csi"):
                Path(str(result.bam) + ext).unlink(missing_ok=True)
        else:
            bam_stats, batch_introns = scan_batch_bam(result.bam, self.config.tile_size)
        spliced_pct = (
            100.0 * bam_stats.n_spliced / bam_stats.n_reads
            if bam_stats.n_reads > 0
            else 0.0
        )

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
        run.obs_version += 1

        run.times_downloaded += br.n_batches
        nd = run.times_downloaded
        w = br.n_batches / nd
        run.avg_umr_pct += (br.uniq_pct - run.avg_umr_pct) * w
        run.avg_spliced_pct += (br.spliced_pct - run.avg_spliced_pct) * w

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
        # The aligner thread may be reading the DB right now (align-ahead):
        # write a temporary file and rename it, so the aligner sees either
        # the old or the new complete file.
        tmp = self._splice_db_path.with_name(self._splice_db_path.name + ".tmp")
        if self.config.longreads:
            n = write_minimap2_junc_bed(introns, tmp)
        else:
            n = write_hisat2_splice_sites(introns, tmp)
        if tmp.exists():
            os.replace(tmp, self._splice_db_path)
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
            threads=self._merge_threads,
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
        self, task: BatchTask, br: BatchResult, t_db: float, t_est: float,
        t_batch: float, t_wait: float = 0.0,
    ) -> None:
        umrs = sum(br.bam_stats.umr_counts.values()) if br.bam_stats else 0
        line = (
            f"{self.batch_count}\t{task.run.record.accession}\t{task.n}\t{task.x}\t"
            f"{int(br.success)}\t{br.uniq_pct:.2f}\t{umrs}\t{task.t_download:.2f}\t"
            f"{br.t_align:.2f}\t{br.t_scan:.2f}\t{t_db:.2f}\t{t_est:.2f}\t"
            f"{t_batch:.2f}\t{len(self._inflight)}\t{int(task.sra_path is not None)}\t"
            f"{t_wait:.2f}\t{task.n_batches}\n"
        )
        log.info(
            "TIMING batch=%d run=%s dl=%.1fs align=%.1fs scan=%.1fs db=%.1fs est=%.1fs "
            "total=%.1fs wait=%.1fs inflight=%d S=%.1f",
            self.batch_count, task.run.record.accession, task.t_download,
            br.t_align, br.t_scan, t_db, t_est, t_batch, t_wait, len(self._inflight),
            self.total_score,
        )
        try:
            with self._timings_path.open("a", encoding="utf-8") as f:
                f.write(line)
        except OSError as e:  # pragma: no cover
            log.debug("could not write timings: %s", e)

    # ------------------------------------------------------------------
    # Output helpers
    # ------------------------------------------------------------------

    def _export_stats(self) -> None:
        """Write per-batch coverage trace and run statistics."""
        stats_path = self.config.outdir / "RunStatistics.csv"

        # Coverage trace: every coverage_trace batches, write a snapshot.
        # A merged download advances batch_count by several, so test whether
        # a multiple was crossed rather than hit exactly.
        prev = getattr(self, "_last_export_count", 0)
        self._last_export_count = self.batch_count
        ct = self.config.coverage_trace
        if ct > 0 and self.batch_count // ct > prev // ct:
            trace = self.config.outdir / f"Coverage{self.batch_count}.tsv"
            write_coverage(self.total_obs, trace)

        # RunStatistics: every batch for the first 10, then every 10
        if self.batch_count < 10 or (self.batch_count - 1) // 10 > (prev - 1) // 10:
            # Pool members are all-zero rows; _finalize writes every run.
            write_run_statistics(
                self._est_runs if self._pool_rep is not None else self.runs, stats_path)

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
    # One seed per run; the batch order itself is built on first pick.
    return [RunState.lazy(r, batch_size, rng.getrandbits(64)) for r in records]


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

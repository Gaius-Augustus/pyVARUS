"""Greedy run selection for ``varus logan``.

Given per-run tile-coverage vectors (from :func:`varus.logan.scan_chunk_bam`)
this module picks the runs that, taken together, cover the genome's
transcribed tiles most evenly. The objective is the same concave "score"
VARUS uses for its sampling controller::

    S(c) = sum_j log(1 + c_j)

where ``c_j`` is the accumulated pseudo-UMR mass on tile ``j``. Every run's
tile vector is normalised to a fixed total mass ``scale`` (roughly the number
of UMRs one VARUS batch of that run would contribute), so a run that spreads
its mass over many tiles has a larger marginal gain than one that piles it on
a handful of rRNA / highly expressed loci. The gain of adding run ``r`` is::

    gain_r = sum_j [ log1p(c_j + l_rj) - log1p(c_j) ]

which is monotone and submodular in the selected set, so the lazy greedy
algorithm of Minoux (1978) gives the same order as naive greedy while
recomputing only a few gains per step.

Only runs with ``status == "accepted"`` participate. The function mutates
``stats[acc].rank`` / ``stats[acc].gain`` in place and returns the coverage
curve so the caller can write it into ``logan_summary.json``.
"""

from __future__ import annotations

import heapq
import logging
from typing import Dict, List, Tuple

import numpy as np

log = logging.getLogger(__name__)


def _sparse_rows(
    stats: Dict[str, "LoganRunStats"],  # noqa: F821 - avoid circular import
    accs: List[str],
    scale: float,
) -> Tuple[Dict[str, Tuple[np.ndarray, np.ndarray]], int]:
    """Return ``{acc: (tile_index_array, mass_array)}`` and the union size.

    Each row is normalised to sum to ``scale``. Rows are kept sparse (index +
    value arrays) instead of a dense R x T matrix; the gain formula only
    depends on non-zero entries, so this is mathematically identical to the
    dense version but avoids a 500 x 200k float32 block for large genomes.
    """
    tile_ids: Dict[Tuple[str, int], int] = {}
    for acc in accs:
        for tile in stats[acc].tiles:
            if tile not in tile_ids:
                tile_ids[tile] = len(tile_ids)
    rows: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for acc in accs:
        tiles = stats[acc].tiles
        if not tiles:
            rows[acc] = (np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.float32))
            continue
        idx = np.fromiter((tile_ids[t] for t in tiles), dtype=np.int64, count=len(tiles))
        val = np.fromiter(tiles.values(), dtype=np.float64, count=len(tiles))
        total = float(val.sum())
        if total > 0:
            # A run only yields reads from the target genome in proportion to
            # its yield (read-mass fraction on target); scale the pseudo-UMRs
            # of one batch accordingly so mixed/contaminated samples with broad
            # but thin coverage do not outrank clean runs.
            y = _yield_factor(stats[acc])
            val = val * (scale * y / total)
        rows[acc] = (idx, val.astype(np.float32))
    return rows, len(tile_ids)


def _yield_factor(st) -> float:
    """Yield fraction in (0, 1]; 1 when the stats carry no yield estimate."""
    y = getattr(st, "yield_pct", None)
    if y is None or y <= 0:
        return 1.0
    return max(0.01, min(1.0, float(y) / 100.0))


def _gain(c: np.ndarray, row: Tuple[np.ndarray, np.ndarray]) -> float:
    idx, val = row
    if idx.size == 0:
        return 0.0
    cj = c[idx]
    return float(np.sum(np.log1p(cj + val) - np.log1p(cj)))


def greedy_select(
    stats: Dict[str, "LoganRunStats"],  # noqa: F821
    *,
    top_k: int,
    scale: float,
    tie_tol: float = 0.01,
    min_gain_frac: float = 1e-3,
) -> List[Tuple[str, float, float, int]]:
    """Lazy-greedy selection of complementary runs.

    Parameters
    ----------
    stats
        ``{acc: LoganRunStats}``; only entries with ``status == "accepted"``
        are candidates.
    top_k
        Maximum number of runs to select (``<= 0`` means no limit).
    scale
        Pseudo-UMR mass of one batch for a run whose reads all map; each
        run's tile vector is scaled to ``scale x yield`` (``yield_pct/100``).
    tie_tol
        Candidates whose fresh gain is within ``tie_tol * best`` of the best
        are considered tied; ties go to a run from a bioproject not selected
        yet, then to the lexicographically smallest accession.
    min_gain_frac
        Stop when the best gain drops below ``min_gain_frac * S`` where ``S``
        is the current objective value. At least one run is always selected
        (if any accepted run exists).

    Returns
    -------
    list of ``(acc, gain, cumulative_S, cumulative_tiles)`` in selection
    order. ``stats[acc].rank`` (1-based) and ``stats[acc].gain`` are set for
    the selected runs; all other accepted runs keep ``rank == 0``.
    """
    accs = sorted(a for a, s in stats.items() if s.status == "accepted")
    for a in accs:
        stats[a].rank = 0
        stats[a].gain = 0.0
    if not accs:
        log.info("greedy_select: no accepted runs")
        return []
    if top_k is None or top_k <= 0:
        top_k = len(accs)

    rows, n_tiles = _sparse_rows(stats, accs, scale)
    c = np.zeros(n_tiles, dtype=np.float64)

    # Max-heap of (-stale_gain, acc). Initially every gain is fresh.
    stale: Dict[str, float] = {a: _gain(c, rows[a]) for a in accs}
    heap: List[Tuple[float, str]] = [(-g, a) for a, g in stale.items()]
    heapq.heapify(heap)
    remaining = set(accs)
    selected_projects: set = set()
    curve: List[Tuple[str, float, float, int]] = []
    S = 0.0

    while remaining and len(curve) < top_k:
        # --- lazy evaluation: find the true maximum ---------------------
        best_acc = None
        best_gain = 0.0
        while heap:
            neg_g, acc = heapq.heappop(heap)
            if acc not in remaining:
                continue
            fresh = _gain(c, rows[acc])
            stale[acc] = fresh
            # Next stale value is an upper bound on every other candidate.
            next_bound = -heap[0][0] if heap else -np.inf
            if fresh >= next_bound:
                best_acc, best_gain = acc, fresh
                break
            heapq.heappush(heap, (-fresh, acc))
        if best_acc is None:
            break

        # --- tie-break among near-equal candidates -----------------------
        thr = (1.0 - tie_tol) * best_gain
        tied: List[Tuple[str, float]] = [(best_acc, best_gain)]
        if tie_tol > 0 and best_gain > 0:
            # Stale gains upper-bound fresh gains (submodularity), so any tied
            # candidate must have a stale value >= thr; refresh only those.
            for a in remaining:
                if a == best_acc or stale[a] < thr:
                    continue
                fresh = _gain(c, rows[a])
                stale[a] = fresh
                if fresh >= thr:
                    tied.append((a, fresh))

        def _key(item: Tuple[str, float]) -> Tuple[int, str]:
            a, _ = item
            bp = stats[a].bioproject or ""
            return (1 if bp in selected_projects else 0, a)

        chosen, chosen_gain = min(tied, key=_key)

        # --- stop rule ---------------------------------------------------
        if curve and chosen_gain < min_gain_frac * S:
            log.info(
                "greedy_select: stopping, best gain %.3g < %.3g x S (%.3g)",
                chosen_gain, min_gain_frac, S,
            )
            break

        # --- accept --------------------------------------------------------
        idx, val = rows[chosen]
        if idx.size:
            c[idx] += val
        S += chosen_gain
        remaining.discard(chosen)
        bp = stats[chosen].bioproject
        if bp:
            selected_projects.add(bp)
        rank = len(curve) + 1
        stats[chosen].rank = rank
        stats[chosen].gain = chosen_gain
        cum_tiles = int(np.count_nonzero(c))
        curve.append((chosen, chosen_gain, S, cum_tiles))
        heap = [(-stale[a], a) for a in remaining]
        heapq.heapify(heap)
        log.info(
            "select %3d %-14s bioproject=%-12s gain=%10.2f S=%12.2f tiles=%d",
            rank, chosen, bp or "-", chosen_gain, S, cum_tiles,
        )
    return curve

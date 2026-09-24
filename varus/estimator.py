"""AdvancedEstimator: tile probability distribution from paper eq. 3.

Stanke et al. (2019) BMC Bioinformatics, eq. 3:

    p̂_r[j] ∝ c^r_j + λ·T·p̄[j] + a

where:
  c^r_j = UMR count from run r in tile j
  p̄[j]  = c_total[j] / Σ c_total  (pooled UMR fraction)
  T      = number of tiles with at least one total observation
  λ      = smoothing coefficient (default 10.0)
  a      = pseudo-count (default 1.0)

Runs not yet downloaded share a common prior that is computed from the pooled
observations of all downloaded runs.  This matches the ``pRep`` optimization
in ``AdvancedEstimator.cpp``: only one p-vector is computed for undownloaded
runs and the rest point to it.

Logan prior (v2 extension)
--------------------------
A run may carry *pseudo-observations* ``ℓ^r_j`` derived from aligning its
Logan contigs to the genome (``varus logan``). They enter eq. 3 exactly like
real counts:

    p̂_r[j] ∝ c^r_j + ℓ^r_j + λ·T·p̄[j] + a

``p̄`` and ``T`` are still computed from real pooled observations only, so a
run with no real batches and no pseudo-observations gets the shared prior as
before, and with ``ℓ = 0`` everywhere the estimator is unchanged.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

Tile = Tuple[str, int]


class AdvancedEstimator:
    """Compute per-run tile probability distributions."""

    def __init__(self, lambda_: float = 10.0, pseudo_count: float = 1.0) -> None:
        self.lambda_ = lambda_
        self.pseudo_count = pseudo_count

    def estimate_sparse(
        self,
        total_arr: np.ndarray,
        run_counts: List[Optional[Tuple[np.ndarray, np.ndarray]]],
    ) -> List[np.ndarray]:
        """Sparse variant: per-run counts as ``(tile indices, values)``.

        ``total_arr`` is the pooled count per tile (length T). Each entry of
        ``run_counts`` holds the run's real plus pseudo counts on the tiles
        it touches (indices into ``total_arr``, unique), or ``None`` for a
        run with neither. Runs with ``None`` all receive the *same* array
        object (the shared prior), so callers can detect sharing with ``is``
        and avoid recomputing per-run quantities for them.

        Cost is O(T) once plus O(T) numpy work per run with counts, instead
        of a Python loop over all T tiles per run.
        """
        T = int(total_arr.shape[0])
        if T == 0:
            return [np.zeros(0) for _ in run_counts]

        total_sum = float(total_arr.sum())
        p_total = total_arr / total_sum if total_sum > 0 else np.full(T, 1.0 / T)
        smooth = self.pseudo_count + self.lambda_ * p_total * T
        prior = smooth / smooth.sum()

        results: List[np.ndarray] = []
        for entry in run_counts:
            if entry is None:
                results.append(prior)
                continue
            idx, val = entry
            raw = smooth.copy()
            raw[idx] += val
            results.append(raw / raw.sum())
        return results

    def estimate_arrays(
        self,
        tiles: List[Tile],
        obs_total: Dict[Tile, int],
        run_obs: List[Dict[Tile, int]],
        times_downloaded: List[int],
        prior_obs: Optional[List[Optional[Dict[Tile, float]]]] = None,
    ) -> List[np.ndarray]:
        """Array variant: one float64 vector per run, aligned to ``tiles``.

        Dict front-end to :meth:`estimate_sparse`; the sharing semantics are
        the same. Tiles of ``run_obs``/``prior_obs`` missing from ``tiles``
        are ignored.
        """
        T = len(tiles)
        if T == 0:
            return [np.zeros(0) for _ in run_obs]
        index = {t: i for i, t in enumerate(tiles)}
        total_arr = np.array(
            [obs_total.get(t, 0) for t in tiles], dtype=np.float64
        )
        entries: List[Optional[Tuple[np.ndarray, np.ndarray]]] = []
        for i, (obs, nd) in enumerate(zip(run_obs, times_downloaded)):
            pseudo = prior_obs[i] if prior_obs is not None else None
            if nd == 0 and not pseudo:
                entries.append(None)
                continue
            entries.append(sparse_counts(index, obs, pseudo))
        return self.estimate_sparse(total_arr, entries)

    def estimate(
        self,
        tiles: List[Tile],
        obs_total: Dict[Tile, int],
        run_obs: List[Dict[Tile, int]],
        times_downloaded: List[int],
        prior_obs: Optional[List[Optional[Dict[Tile, float]]]] = None,
    ) -> List[Dict[Tile, float]]:
        """Return one {tile: probability} dict per run.

        Parameters
        ----------
        tiles:            Ordered list of all tiles with nonzero pooled count
                          (plus any tiles that only appear in ``prior_obs``).
        obs_total:        Pooled UMR counts across every run.
        run_obs:          Per-run {tile: count} observations.
        times_downloaded: Download count per run (parallel to run_obs).
        prior_obs:        Optional per-run {tile: pseudo-count} (Logan prior).
        """
        if not tiles:
            return [{} for _ in run_obs]
        arrays = self.estimate_arrays(
            tiles, obs_total, run_obs, times_downloaded, prior_obs
        )
        results: List[Dict[Tile, float]] = []
        shared: Optional[Dict[Tile, float]] = None
        shared_arr = None
        for arr in arrays:
            # Share one dict for all runs on the shared prior (pRep semantics).
            if shared_arr is not None and arr is shared_arr:
                results.append(shared)  # type: ignore[arg-type]
                continue
            d = dict(zip(tiles, arr.tolist()))
            if shared_arr is None and _is_shared_candidate(arr, arrays):
                shared_arr, shared = arr, d
            results.append(d)
        return results


def sparse_counts(
    index: Dict[Tile, int],
    obs: Optional[Dict[Tile, int]],
    pseudo: Optional[Dict[Tile, float]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Merge real and pseudo counts of one run into ``(indices, values)``.

    Indices are unique and refer to ``index``; tiles not in ``index`` are
    dropped. Cost is proportional to the run's own tiles, not to T.
    """
    merged: Dict[int, float] = {}
    for src in (obs, pseudo):
        if not src:
            continue
        for t, v in src.items():
            i = index.get(t)
            if i is not None:
                merged[i] = merged.get(i, 0.0) + float(v)
    if not merged:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.float64)
    idx = np.fromiter(merged.keys(), dtype=np.int64, count=len(merged))
    val = np.fromiter(merged.values(), dtype=np.float64, count=len(merged))
    return idx, val


def _is_shared_candidate(arr: np.ndarray, arrays: List[np.ndarray]) -> bool:
    """True when ``arr`` is the object handed to more than one run."""
    n = 0
    for a in arrays:
        if a is arr:
            n += 1
            if n > 1:
                return True
    return False

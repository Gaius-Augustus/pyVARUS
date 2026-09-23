"""Tests for varus.select.greedy_select (toy tile matrices, no I/O)."""

from __future__ import annotations

import math
import random
from typing import Dict, List

import numpy as np
import pytest

from varus.logan import LoganRunStats
from varus.select import greedy_select


def _stats_from_rows(rows: Dict[str, List[float]], bioprojects=None,
                     status="accepted") -> Dict[str, LoganRunStats]:
    stats = {}
    for acc, row in rows.items():
        tiles = {("chr1", j): float(v) for j, v in enumerate(row) if v > 0}
        stats[acc] = LoganRunStats(acc=acc, n_contigs=1000, n_tiles=len(tiles),
                                   tile_mass=sum(tiles.values()), tiles=tiles,
                                   status=status,
                                   bioproject=(bioprojects or {}).get(acc, ""))
    return stats


def _brute_force(rows: Dict[str, List[float]], scale: float, top_k: int,
                 min_gain_frac: float = 1e-3) -> List[str]:
    """Naive greedy with full recomputation each step (dense matrix)."""
    accs = sorted(rows)
    T = len(next(iter(rows.values())))
    L = np.zeros((len(accs), T))
    for i, a in enumerate(accs):
        r = np.asarray(rows[a], dtype=float)
        L[i] = r * (scale / r.sum()) if r.sum() > 0 else r
    c = np.zeros(T)
    chosen: List[str] = []
    S = 0.0
    remaining = list(range(len(accs)))
    while remaining and len(chosen) < top_k:
        gains = [(np.sum(np.log1p(c + L[i]) - np.log1p(c)), accs[i], i) for i in remaining]
        best = max(t[0] for t in gains)
        # exact ties: smallest acc
        tied = sorted(t for t in gains if math.isclose(t[0], best, rel_tol=1e-9))
        g, a, i = tied[0]
        if chosen and g < min_gain_frac * S:
            break
        c += L[i]
        S += g
        chosen.append(a)
        remaining.remove(i)
    return chosen


def test_greedy_matches_brute_force_on_random_matrix():
    rng = random.Random(7)
    rows: Dict[str, List[float]] = {}
    for i in range(12):
        # sparse random rows with varying breadth
        row = [rng.choice([0, 0, 0, 1, 3, 10]) * rng.random() for _ in range(40)]
        if sum(row) == 0:
            row[i] = 1.0
        rows[f"R{i:02d}"] = row
    stats = _stats_from_rows(rows)
    curve = greedy_select(stats, top_k=8, scale=100.0, tie_tol=0.0)
    got = [acc for acc, *_ in curve]
    assert got == _brute_force(rows, 100.0, 8)
    # ranks set in place, 1-based, contiguous
    assert [stats[a].rank for a in got] == list(range(1, len(got) + 1))
    for a in stats:
        if a not in got:
            assert stats[a].rank == 0
    # cumulative_S is increasing and gains are positive
    S = [s for _, _, s, _ in curve]
    assert all(b > a for a, b in zip(S, S[1:]))
    assert all(g > 0 for _, g, _, _ in curve)
    assert math.isclose(S[-1], sum(g for _, g, _, _ in curve), rel_tol=1e-6)


def test_broad_run_beats_concentrated_run():
    rows = {
        "BROAD": [1.0] * 20,
        "PEAK": [100.0] + [0.0] * 19,
    }
    stats = _stats_from_rows(rows)
    curve = greedy_select(stats, top_k=2, scale=50.0)
    assert curve[0][0] == "BROAD"
    assert curve[0][3] == 20  # cumulative tiles


def test_complementary_run_preferred_over_duplicate():
    rows = {
        "A": [1.0] * 10 + [0.0] * 10,
        "A2": [1.0] * 10 + [0.0] * 10,  # identical to A
        "B": [0.0] * 10 + [1.0] * 10,   # complementary
    }
    stats = _stats_from_rows(rows)
    curve = greedy_select(stats, top_k=3, scale=50.0)
    order = [acc for acc, *_ in curve]
    assert order[0] == "A"          # tie between A and A2 -> lexicographic
    assert order[1] == "B"          # complementary before duplicate
    assert order[2] == "A2"


def test_tie_break_prefers_unselected_bioproject():
    rows = {
        "SRR1": [1.0] * 10,
        "SRR2": [1.0] * 10,   # identical gain to SRR1 and SRR3
        "SRR3": [1.0] * 10,
    }
    bps = {"SRR1": "PRJA", "SRR2": "PRJA", "SRR3": "PRJB"}
    stats = _stats_from_rows(rows, bps)
    curve = greedy_select(stats, top_k=3, scale=10.0)
    order = [acc for acc, *_ in curve]
    # first pick: all tied, no project selected yet -> lexicographic SRR1
    assert order[0] == "SRR1"
    # second pick: SRR2 (PRJA, already selected) vs SRR3 (PRJB, new) -> SRR3
    assert order[1] == "SRR3"
    assert order[2] == "SRR2"


def test_tie_tol_near_equal_gains():
    """A slightly worse run from a new bioproject wins within tie_tol."""
    rows = {
        "SRR0": [5.0] * 50 + [0.0] * 51,       # clearly best first pick
        "SRR1": [1.0] * 100 + [0.0],           # PRJX (same as SRR0)
        "SRR2": [1.0] * 99 + [0.0, 1.0],       # PRJY, same breadth as SRR1
    }
    rows["SRR1"][0] = 1.02  # SRR1 is a hair better than SRR2 on its own
    bps = {"SRR0": "PRJX", "SRR1": "PRJX", "SRR2": "PRJY"}
    stats = _stats_from_rows(rows, bps)
    curve = greedy_select(stats, top_k=2, scale=10.0, tie_tol=0.05)
    assert [a for a, *_ in curve] == ["SRR0", "SRR2"]
    # with a zero tolerance the strictly better run wins regardless of project
    stats = _stats_from_rows(rows, bps)
    curve = greedy_select(stats, top_k=2, scale=10.0, tie_tol=0.0)
    assert [a for a, *_ in curve][1] in ("SRR1", "SRR2")
    g1 = _brute_force(rows, 10.0, 2)
    assert [a for a, *_ in curve] == g1


def test_top_k_limits_selection():
    rows = {f"R{i}": [1.0 if j % 5 == i else 0.0 for j in range(25)] for i in range(5)}
    stats = _stats_from_rows(rows)
    curve = greedy_select(stats, top_k=2, scale=10.0)
    assert len(curve) == 2
    assert sum(1 for s in stats.values() if s.rank) == 2
    curve = greedy_select(stats, top_k=0, scale=10.0)  # 0 = unlimited
    assert len(curve) == 5


def test_stop_rule_min_gain_frac():
    """A duplicate of an already-selected run adds little; stop when the
    marginal gain falls below ``min_gain_frac * S``."""
    rows = {
        "BIG": [1.0] * 1000,
        "DUP": [1.0] * 1000,
    }
    stats = _stats_from_rows(rows)
    # scale 1e6 -> c_j = 1000 after BIG; S = 1000*log1p(1000) ~ 6908,
    # DUP gain = 1000*(log1p(2000) - log1p(1000)) ~ 693 (~10 % of S)
    curve = greedy_select(stats, top_k=5, scale=1e6, min_gain_frac=0.2)
    assert [a for a, *_ in curve] == ["BIG"]
    assert stats["DUP"].rank == 0
    assert curve[0][2] == pytest.approx(1000 * math.log1p(1000), rel=1e-6)
    # with a lower threshold both get picked
    curve = greedy_select(stats, top_k=5, scale=1e6, min_gain_frac=0.05)
    assert [a for a, *_ in curve] == ["BIG", "DUP"]
    assert curve[1][1] == pytest.approx(1000 * (math.log1p(2000) - math.log1p(1000)), rel=1e-6)


def test_always_selects_at_least_one_and_ignores_non_accepted():
    rows = {"A": [1.0] * 5, "B": [1.0] * 50}
    stats = _stats_from_rows(rows)
    stats["B"].status = "rejected"
    curve = greedy_select(stats, top_k=10, scale=1.0)
    assert [a for a, *_ in curve] == ["A"]
    assert stats["B"].rank == 0


def test_no_accepted_runs_returns_empty():
    stats = _stats_from_rows({"A": [1.0]}, status="rejected")
    assert greedy_select(stats, top_k=3, scale=1.0) == []


def test_run_without_tiles_gets_zero_gain_and_is_not_first():
    rows = {"EMPTY": [0.0] * 10, "OK": [1.0] * 10}
    stats = _stats_from_rows(rows)
    curve = greedy_select(stats, top_k=2, scale=10.0, min_gain_frac=0.0)
    assert curve[0][0] == "OK"
    assert curve[1][0] == "EMPTY" and curve[1][1] == 0.0

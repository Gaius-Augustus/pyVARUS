"""Tests for the v2 speed-ups: incremental strand DB, aligner flags, rolling
merge, prefetch, parallel downloads with in-flight accounting, exit status,
Logan prior wiring. Everything is mocked; no network, no aligners."""

from __future__ import annotations

import math
import random
import threading
import subprocess
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from varus import align, download, merge
from varus.controller import (
    EXIT_NO_USABLE_DATA, BatchTask, Controller, RunState, VARUSConfig,
    apply_logan_prior,
)
from varus.download import BatchPaths
from varus.estimator import AdvancedEstimator
from varus.introns import IntronCounts
from varus.runlist import RunRecord
from varus.tiles import BAMStats

try:
    from tests.conftest import requires_pysam
except ImportError:  # pragma: no cover
    requires_pysam = pytest.mark.skip(reason="conftest not found")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _rec(acc="SRR1", spots=1_000_000, avg_len=100.0, bioproject="") -> RunRecord:
    return RunRecord(accession=acc, total_spots=spots, total_bases=int(spots * avg_len),
                     avg_len=avg_len, paired=False, colorspace=False,
                     platform="ILLUMINA", bioproject=bioproject)


def _cfg(tmp_path: Path, **kw) -> VARUSConfig:
    d = dict(genome=tmp_path / "g.fa", index_prefix=tmp_path / "idx",
             outdir=tmp_path / "out", batch_size=50_000, max_batches=6,
             tile_size=5_000, merge_every=0, hisat2_mm=True,
             parallel_downloads=1)  # serial unless a test asks otherwise
    d.update(kw)
    return VARUSConfig(**d)


def _write_genome(path: Path) -> None:
    # chr1: 'gt' at 10..11 (donor of intron starting at 11, 1-based) and 'ag'
    # ending at 50; chr2 has a non-canonical site.
    seq = list("a" * 100)
    seq[10:12] = "gt"      # 0-based 10,11 -> 1-based 11,12 = intron start 11
    seq[48:50] = "ag"      # 0-based 48,49 -> 1-based 49,50 = intron end 50
    path.write_text(">chr1\n" + "".join(seq) + "\n>chr2\n" + "c" * 100 + "\n")


class _FakePipe:
    def close(self):
        pass


class _FakeProc:
    def __init__(self, cmd, **kw):
        self.cmd = list(cmd)
        self.stdout = _FakePipe()
        self.returncode = 0

    def wait(self):
        return 0


def _mock_stack(monkeypatch, tile_by_run, uniq=50.0, introns=None):
    """Patch download/align/scan so Controller.run() works without tools.

    ``tile_by_run``: acc -> dict of tile->count returned for every batch of
    that run. A run with uniq below the gate can be simulated by passing
    ``uniq`` as a dict acc->uniq_pct.
    """
    calls = {"download": [], "align": []}

    def fake_download(accession, n, x, paired, outdir, sra_path=None, **kw):
        bdir = download.batch_dir_for(outdir, accession, n, x)
        bdir.mkdir(parents=True, exist_ok=True)
        r1 = bdir / f"{accession}.fasta"
        r1.write_text(">r\nACGT\n")
        calls["download"].append((accession, n, x, sra_path))
        return BatchPaths(r1=r1, r2=None, batch_dir=bdir)

    def fake_align(r1, r2, *, index_prefix, batch_dir, threads, intron_db=None, **kw):
        bam = batch_dir / "Aligned.out.bam"
        bam.write_bytes(b"BAM")
        logp = batch_dir / "Log.final.out"
        logp.write_text("x")
        calls["align"].append((batch_dir, intron_db, kw))
        return align.AlignmentResult(bam=bam, log=logp)

    def fake_parse(log_path, batch_size):
        acc = log_path.parent.parent.name
        u = uniq.get(acc, 50.0) if isinstance(uniq, dict) else uniq
        return {"num_uniq": u * batch_size / 100, "uniq_pct": u}

    def fake_scan(bam, tile_size):
        acc = bam.parent.parent.name
        counts = dict(tile_by_run[acc])
        stats = BAMStats(umr_counts=counts, n_reads=sum(counts.values()) or 1,
                         n_spliced=0)
        return stats, IntronCounts(dict(introns or {}))

    monkeypatch.setattr("varus.controller.download_batch", fake_download)
    monkeypatch.setattr("varus.controller.align_batch_hisat2", fake_align)
    monkeypatch.setattr("varus.controller.parse_hisat2_log", fake_parse)
    monkeypatch.setattr("varus.controller.scan_batch_bam", fake_scan)
    monkeypatch.setattr("varus.controller.scan_batch_bam_parallel",
                        lambda bam, tile_size, ex, n_parts: fake_scan(bam, tile_size))
    merged = []

    def fake_merge(bams, out, threads=4, compression=None, **kw):
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(b"MERGED")
        merged.append((list(bams), out, compression))
        return out

    monkeypatch.setattr("varus.controller.merge_bams", fake_merge)
    calls["merged"] = merged
    return calls


# ---------------------------------------------------------------------------
# StrandAssigner
# ---------------------------------------------------------------------------

@requires_pysam
def test_strand_assigner_incremental(tmp_path: Path):
    from varus.strand import StrandAssigner, assign_strand
    g = tmp_path / "g.fa"
    _write_genome(g)
    sa = StrandAssigner(g)
    batch = IntronCounts({("chr1", 11, 50, "."): 3, ("chr2", 11, 50, "."): 1})
    stranded, n_new = sa.assign_new(batch)
    assert n_new == 2
    assert stranded.counts == {("chr1", 11, 50, "+"): 3}   # chr2 non-canonical dropped
    assert sa.n_cached == 2
    # Second batch: same keys are cache hits, one new key.
    batch2 = IntronCounts({("chr1", 11, 50, "."): 2, ("chr1", 11, 40, "."): 1})
    stranded2, n_new2 = sa.assign_new(batch2)
    assert n_new2 == 1
    assert stranded2.counts[("chr1", 11, 50, "+")] == 2
    # The functional wrapper still works and agrees.
    assert assign_strand(batch, g).counts == stranded.counts
    # preload seeds the cache without touching the genome.
    sa2 = StrandAssigner(tmp_path / "missing.fa")
    assert sa2.preload(IntronCounts({("chrX", 5, 9, "-"): 1, ("chrX", 5, 9, "."): 1})) == 1
    assert sa2.resolve("chrX", 5, 9) == "-"


def test_rebuild_intron_db_only_rewrites_on_new_keys(tmp_path: Path):
    cfg = _cfg(tmp_path)
    ctrl = Controller(cfg, [])
    ctrl.stranded_introns = IntronCounts({("chr1", 10, 50, "+"): 1})
    ctrl.cumulative_introns = IntronCounts({("chr1", 10, 50, "."): 1})
    ctrl._db_new_keys = 1
    with patch("varus.controller.write_hisat2_splice_sites", return_value=1) as m:
        ctrl._rebuild_intron_db()
        assert m.call_count == 1
        ctrl._rebuild_intron_db()          # nothing new -> no rewrite
        assert m.call_count == 1
        ctrl._db_new_keys = 2
        ctrl._rebuild_intron_db()
        assert m.call_count == 2


def test_rebuild_intron_db_min_mult_filter(tmp_path: Path):
    cfg = _cfg(tmp_path, splice_db_min_mult=2)
    ctrl = Controller(cfg, [])
    ctrl.stranded_introns = IntronCounts({("chr1", 10, 50, "+"): 1, ("chr1", 60, 90, "-"): 3})
    ctrl.cumulative_introns = IntronCounts({("chr1", 10, 50, "."): 1})
    ctrl._db_new_keys = 2
    with patch("varus.controller.write_hisat2_splice_sites", return_value=1) as m:
        ctrl._rebuild_intron_db()
    written = m.call_args.args[0]
    assert list(written.counts) == [("chr1", 60, 90, "-")]


# ---------------------------------------------------------------------------
# align / merge / download flags
# ---------------------------------------------------------------------------

def test_hisat2_flags_default_and_off(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(align.shutil, "which", lambda t: f"/fake/{t}")
    cmds = []
    monkeypatch.setattr(align.subprocess, "Popen",
                        lambda cmd, **kw: cmds.append(list(cmd)) or _FakeProc(cmd, **kw))
    align.align_batch_hisat2(r1=tmp_path / "r.fa", r2=None, index_prefix=tmp_path / "i",
                             batch_dir=tmp_path / "b", threads=4)
    h, s = cmds[0], cmds[1]
    assert "--mm" in h and "--no-unal" in h
    assert s[s.index("-l") + 1] == "1"
    cmds.clear()
    align.align_batch_hisat2(r1=tmp_path / "r.fa", r2=None, index_prefix=tmp_path / "i",
                             batch_dir=tmp_path / "b", threads=4, mm=False,
                             keep_unaligned=True, sort_compression=None)
    h, s = cmds[0], cmds[1]
    assert "--mm" not in h and "--no-unal" not in h and "-l" not in s


def test_align_contigs_minimap2_command(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(align.shutil, "which", lambda t: f"/fake/{t}")
    cmds = []
    monkeypatch.setattr(align.subprocess, "Popen",
                        lambda cmd, **kw: cmds.append(list(cmd)) or _FakeProc(cmd, **kw))
    out = align.align_contigs_minimap2([tmp_path / "a.fa", tmp_path / "b.fa"],
                                       index=tmp_path / "g.mmi", out_bam=tmp_path / "o/x.bam",
                                       threads=3, max_intron=1234)
    mm2 = cmds[0]
    assert mm2[0].endswith("minimap2") and mm2[1:3] == ["-t", "3"]
    assert "--secondary=no" in mm2 and mm2[mm2.index("-G") + 1] == "1234"
    assert mm2[-2:] == [str(tmp_path / "a.fa"), str(tmp_path / "b.fa")]
    assert out == tmp_path / "o/x.bam"
    with pytest.raises(ValueError):
        align.align_contigs_minimap2([], index=tmp_path / "g.mmi", out_bam=tmp_path / "y.bam")


def test_merge_bams_compression_flag(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(merge.shutil, "which", lambda t: "/fake/samtools")
    captured = {}
    monkeypatch.setattr(merge.subprocess, "run",
                        lambda cmd, check: captured.setdefault("cmd", list(cmd)))
    merge.merge_bams([tmp_path / "a.bam"], tmp_path / "o.bam", threads=2, compression=1)
    cmd = captured["cmd"]
    assert cmd[cmd.index("-l") + 1] == "1" and cmd[cmd.index("-@") + 1] == "2"


def test_download_batch_uses_local_sra(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(download.shutil, "which", lambda _: "/fake/fastq-dump")
    captured = {}

    def fake_run(cmd, check):
        captured["cmd"] = list(cmd)
        bdir = Path(cmd[cmd.index("-O") + 1])
        (bdir / "SRR9.fasta").write_text("")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(download.subprocess, "run", fake_run)
    sra = tmp_path / "sra" / "SRR9" / "SRR9.sra"
    paths = download.download_batch("SRR9", 0, 9, paired=False, outdir=tmp_path, sra_path=sra)
    assert captured["cmd"][-1] == str(sra)
    assert paths.r1.name == "SRR9.fasta"


def test_prefetch_run_command_and_reuse(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(download.shutil, "which", lambda _: "/fake/prefetch")
    captured = []

    def fake_run(cmd, check):
        captured.append(list(cmd))
        acc = cmd[-1]
        d = Path(cmd[cmd.index("-O") + 1]) / acc
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{acc}.sra").write_bytes(b"x")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(download.subprocess, "run", fake_run)
    p = download.prefetch_run("SRR5", tmp_path / "sra", max_size_gb=12.7)
    assert p == tmp_path / "sra" / "SRR5" / "SRR5.sra"
    assert captured[0][captured[0].index("--max-size") + 1] == "12G"
    # second call reuses the file, no subprocess
    assert download.prefetch_run("SRR5", tmp_path / "sra") == p
    assert len(captured) == 1
    assert download.find_prefetched("NOPE", tmp_path / "sra") is None


# ---------------------------------------------------------------------------
# estimator prior
# ---------------------------------------------------------------------------

def test_estimator_prior_shifts_mass_and_degrades():
    est = AdvancedEstimator(lambda_=10.0, pseudo_count=1.0)
    tiles = [("c", 0), ("c", 1), ("c", 2)]
    total = {("c", 0): 10, ("c", 1): 10}
    run_obs = [{}, {}, {("c", 0): 20}]
    nd = [0, 0, 1]
    base = est.estimate(tiles, total, run_obs, nd)
    withp = est.estimate(tiles, total, run_obs, nd,
                         prior_obs=[None, {("c", 2): 1000.0}, None])
    # run 0: unchanged shared prior; run 1: mass moved to tile 2; run 2 unchanged
    assert base[0] == withp[0]
    assert withp[1][("c", 2)] > 0.9 and base[1][("c", 2)] < 0.2
    assert base[2] == withp[2]
    # zero prior == no prior
    zero = est.estimate(tiles, total, run_obs, nd, prior_obs=[None, {}, None])
    assert zero == base
    # Prior washes out as real counts grow.
    strong = est.estimate(tiles, total, [{}, {("c", 0): 10**6}, {}], [0, 1, 0],
                          prior_obs=[None, {("c", 2): 1000.0}, None])
    assert strong[1][("c", 0)] > 0.99


def test_estimator_arrays_share_prior_object():
    est = AdvancedEstimator()
    tiles = [("c", 0)]
    arrs = est.estimate_arrays(tiles, {("c", 0): 1}, [{}, {}, {("c", 0): 1}], [0, 0, 1])
    assert arrs[0] is arrs[1] and arrs[2] is not arrs[0]


# ---------------------------------------------------------------------------
# Controller: profit with priors, lazy greedy, in-flight accounting
# ---------------------------------------------------------------------------

def test_calculate_profit_distinguishes_prior_runs(tmp_path: Path):
    cfg = _cfg(tmp_path)
    rng = random.Random(0)
    a = RunState.from_record(_rec("A"), cfg.batch_size, rng)
    b = RunState.from_record(_rec("B"), cfg.batch_size, rng)
    c = RunState.from_record(_rec("C"), cfg.batch_size, rng)
    a.prior_obs = {("chr1", 0): 100.0}          # tile already well covered
    b.prior_obs = {("chr1", 9): 100.0}          # fresh tile
    ctrl = Controller(cfg, [a, b, c])
    ctrl.total_obs = {("chr1", 0): 500}
    ctrl._estimate_p()
    ctrl._calculate_profit()
    assert b.expected_profit > a.expected_profit
    assert c.expected_profit != b.expected_profit
    # Without priors, undownloaded runs share one profit (legacy shortcut).
    a.prior_obs = {}
    b.prior_obs = {}
    ctrl._logan_tiles = set()
    ctrl._estimate_p()
    ctrl._calculate_profit()
    assert a.expected_profit == b.expected_profit == c.expected_profit


def test_logan_prior_first_only(tmp_path: Path):
    """With logan_prior_first_only a downloaded run's prior is dropped."""
    for first_only in (False, True):
        cfg = _cfg(tmp_path)
        cfg.logan_prior_first_only = first_only
        rng = random.Random(0)
        a = RunState.from_record(_rec("A"), cfg.batch_size, rng)
        b = RunState.from_record(_rec("B"), cfg.batch_size, rng)
        a.prior_obs = {("chr1", 9): 1000.0}
        b.prior_obs = {("chr1", 9): 1000.0}
        ctrl = Controller(cfg, [a, b])
        ctrl.total_obs = {("chr1", 0): 500, ("chr1", 9): 1}
        a.times_downloaded = 1
        a.observations = {("chr1", 0): 500}
        a.obs_version += 1
        ctrl._estimate_p()
        i9 = ctrl._tile_index[("chr1", 9)]
        # b (never downloaded) keeps its prior either way
        assert b.p_arr[i9] > 0.5
        if first_only:
            assert a.p_arr[i9] < 0.01
        else:
            assert a.p_arr[i9] > 0.3


def test_lazy_greedy_matches_brute_force(tmp_path: Path):
    cfg = _cfg(tmp_path)
    rng = random.Random(3)
    runs = [RunState.from_record(_rec(f"R{i}"), cfg.batch_size, rng) for i in range(6)]
    for i, r in enumerate(runs):
        r.prior_obs = {("chr1", j): float((i * 7 + j * 3) % 11 + 1) for j in range(8)}
    ctrl = Controller(cfg, runs)
    ctrl.total_obs = {("chr1", j): 40 * (j % 3) for j in range(8)}
    ctrl._estimate_p()
    ctrl._calculate_profit()
    ctrl._sim_extra = {("chr1", 1): 3000.0, ("chr1", 5): 500.0}
    picked = ctrl._choose_next_run_simulated()
    brute = max(ctrl.downloadable, key=lambda r: ctrl._profit(r, ctrl._sim_extra))
    assert picked is brute
    assert abs(ctrl.max_profit - ctrl._profit(brute, ctrl._sim_extra)) < 1e-9
    # Without anything in flight the simulated pick is the plain pick.
    ctrl._sim_extra = {}
    assert ctrl._choose_next_run_simulated() is ctrl._choose_next_run()


def test_refill_keeps_k_in_flight_and_accounts_expected(tmp_path: Path, monkeypatch):
    cfg = _cfg(tmp_path, parallel_downloads=3, max_batches=10)
    rng = random.Random(1)
    runs = [RunState.from_record(_rec(f"R{i}"), cfg.batch_size, rng) for i in range(4)]
    tiles = {f"R{i}": {("chr1", i): 100, ("chr1", 4): 10} for i in range(4)}
    calls = _mock_stack(monkeypatch, tiles)
    ctrl = Controller(cfg, runs)
    # Give the estimator something to work with, as after one applied batch.
    ctrl.total_obs = {("chr1", 4): 10}
    for r in runs:
        r.prior_obs = {t: float(v) for t, v in tiles[r.record.accession].items()}
    ctrl._logan_tiles = {t for d in tiles.values() for t in d}
    ctrl._estimate_p()
    ctrl._calculate_profit()
    from concurrent.futures import ThreadPoolExecutor
    ctrl._serial = False
    ctrl._max_inflight = 3
    ctrl._dl_ex = ThreadPoolExecutor(max_workers=3)
    try:
        ctrl._refill()
        assert len(ctrl._inflight) == 3
        accs = {t.run.record.accession for t in ctrl._inflight}
        assert len(accs) == 3, "parallel picks must not blindly repeat one run"
        expected_sum = {}
        for t in ctrl._inflight:
            for k, v in t.expected.items():
                expected_sum[k] = expected_sum.get(k, 0.0) + v
        for k, v in expected_sum.items():
            assert abs(ctrl._sim_extra[k] - v) < 1e-9
        task = ctrl._next_ready_task()
        ctrl._inflight.remove(task)
        ctrl._release_expected(task)
        for k, v in task.expected.items():
            assert abs(ctrl._sim_extra.get(k, 0.0) - (expected_sum[k] - v)) < 1e-9
    finally:
        ctrl._drain_inflight()
        ctrl._dl_ex.shutdown(wait=True)


# ---------------------------------------------------------------------------
# Controller.run(): serial determinism, exit status, cleanup, rolling merge
# ---------------------------------------------------------------------------

def _run_ctrl(tmp_path, monkeypatch, tiles, **cfgkw):
    cfg = _cfg(tmp_path, **cfgkw)
    rng = random.Random(cfg.seed if cfg.seed is not None else 0)
    runs = [RunState.from_record(_rec(acc), cfg.batch_size, rng) for acc in tiles]
    calls = _mock_stack(monkeypatch, tiles)
    ctrl = Controller(cfg, runs)
    status = ctrl.run()
    return ctrl, status, calls


def test_run_serial_is_deterministic_and_greedy(tmp_path: Path, monkeypatch):
    tiles = {"A": {("chr1", 0): 100}, "B": {("chr1", 1): 100}, "C": {("chr1", 0): 100}}
    ctrl1, st1, c1 = _run_ctrl(tmp_path / "1", monkeypatch, tiles, seed=7, max_batches=6)
    ctrl2, st2, c2 = _run_ctrl(tmp_path / "2", monkeypatch, tiles, seed=7, max_batches=6)
    assert st1 == st2 == 0
    picks1 = [a for a, *_ in c1["download"]]
    picks2 = [a for a, *_ in c2["download"]]
    assert picks1 == picks2 and len(picks1) == 6
    # B covers a tile nobody else does: it must be picked at least once early
    assert "B" in picks1[:3]
    assert (tmp_path / "1" / "out" / "VARUS.bam").is_file()
    assert (tmp_path / "1" / "out" / "BatchTimings.tsv").read_text().count("\n") == 7
    assert ctrl1.batch_count == 6


def test_run_parallel_completes_and_uses_all_runs(tmp_path: Path, monkeypatch):
    tiles = {f"R{i}": {("chr1", i): 100} for i in range(5)}
    ctrl, st, calls = _run_ctrl(tmp_path, monkeypatch, tiles, seed=1, max_batches=8,
                                parallel_downloads=3)
    assert st == 0 and ctrl.batch_count == 8
    assert len({a for a, *_ in calls["download"]}) >= 4
    assert not ctrl._inflight


def test_run_exit_status_when_all_rejected(tmp_path: Path, monkeypatch):
    cfg = _cfg(tmp_path, max_batches=5)
    rng = random.Random(0)
    runs = [RunState.from_record(_rec(a), cfg.batch_size, rng) for a in ("X", "Y")]
    _mock_stack(monkeypatch, {"X": {}, "Y": {}}, uniq=1.0)
    ctrl = Controller(cfg, runs)
    status = ctrl.run()
    assert status == EXIT_NO_USABLE_DATA
    assert not (cfg.outdir / "VARUS.bam").exists()
    assert (cfg.outdir / "RunStatistics.csv").is_file()
    # Rejected batch directories are removed entirely (no stray BAM/log).
    assert not (cfg.outdir / "batches").exists()


def test_rolling_merge_parts_and_final(tmp_path: Path, monkeypatch):
    tiles = {"A": {("chr1", 0): 100}, "B": {("chr1", 1): 100}}
    ctrl, st, calls = _run_ctrl(tmp_path, monkeypatch, tiles, seed=0, max_batches=5,
                                merge_every=2)
    assert st == 0
    parts = [m for m in calls["merged"] if "part_" in m[1].name]
    final = [m for m in calls["merged"] if m[1].name == "VARUS.bam"]
    assert len(parts) == 2 and all(len(p[0]) == 2 and p[2] == 1 for p in parts)
    assert len(final) == 1
    # final merge takes the 2 parts plus the 1 remaining batch BAM
    assert len(final[0][0]) == 3 and final[0][2] is None
    assert not (tmp_path / "out" / "batches").exists()
    assert not (tmp_path / "out" / "merged").exists()


def test_prefetch_trigger_passes_local_sra(tmp_path: Path, monkeypatch):
    tiles = {"A": {("chr1", 0): 100}, "B": {("chr1", 1): 100}}
    sra = tmp_path / "out" / "sra" / "A" / "A.sra"

    def fake_prefetch(acc, sra_dir, *, max_size_gb, **kw):
        p = sra_dir / acc / f"{acc}.sra"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"sra")
        return p

    monkeypatch.setattr("varus.controller.prefetch_run", fake_prefetch)
    ctrl, st, calls = _run_ctrl(tmp_path, monkeypatch, tiles, seed=0, max_batches=8,
                                prefetch=True, prefetch_after=2, parallel_downloads=1)
    assert st == 0
    local = [(a, s) for a, n, x, s in calls["download"] if s is not None]
    assert local, "later batches of a repeatedly picked run must use the local .sra"
    assert all(s.name == f"{a}.sra" for a, s in local)
    assert not (tmp_path / "out" / "sra").exists()      # scratch removed at the end


def test_prefetch_respects_size_cap(tmp_path: Path):
    cfg = _cfg(tmp_path, prefetch=True, prefetch_after=1, prefetch_max_gb=0.001)
    rng = random.Random(0)
    big = RunState.from_record(_rec("BIG", spots=50_000_000, avg_len=150.0), cfg.batch_size, rng)
    ctrl = Controller(cfg, [big])
    from concurrent.futures import ThreadPoolExecutor
    ctrl._pf_ex = ThreadPoolExecutor(max_workers=1)
    try:
        big.times_downloaded = 1
        ctrl._maybe_prefetch(big)
        assert big.prefetch_future is None and big.prefetch_failed
    finally:
        ctrl._pf_ex.shutdown(wait=True)


# ---------------------------------------------------------------------------
# Logan prior wiring
# ---------------------------------------------------------------------------

def test_apply_logan_prior_filters_and_normalises(tmp_path: Path):
    cfg = _cfg(tmp_path)
    rng = random.Random(0)
    runs = [RunState.from_record(_rec(a), cfg.batch_size, rng) for a in ("A", "B", "C", "D")]
    logan = SimpleNamespace(
        status={"A": "accepted", "B": "accepted", "C": "rejected"},
        rank={"A": 1, "B": 2},
        tiles={"A": {("chr1", 0): 3.0, ("chr1", 1): 1.0}, "B": {("chr1", 5): 2.0}},
        introns=None, splice_sites=None, junc_bed=None, params={},
    )
    kept = apply_logan_prior(runs, logan, batch_size=50_000, prior_batches=1.0)
    assert [r.record.accession for r in kept] == ["A", "B", "D"]
    a = kept[0]
    assert abs(sum(a.prior_obs.values()) - 25_000) < 1e-6
    assert abs(a.prior_obs[("chr1", 0)] / a.prior_obs[("chr1", 1)] - 3.0) < 1e-9
    assert kept[2].prior_obs == {}
    top = apply_logan_prior(
        [RunState.from_record(_rec(a), cfg.batch_size, rng) for a in ("A", "B", "C", "D")],
        logan, batch_size=50_000, prior_batches=0.0, top=1, only=True)
    assert [r.record.accession for r in top] == ["A"] and top[0].prior_obs == {}


def test_controller_seeds_from_logan(tmp_path: Path):
    from varus.introns import write_introns_gff
    ld = tmp_path / "logan"
    ld.mkdir()
    gff = ld / "logan_introns.gff"
    write_introns_gff(IntronCounts({("chr1", 10, 50, "+"): 4}), gff)
    ss = ld / "logan.splice_sites"
    ss.write_text("chr1\t8\t50\t+\n")
    logan = SimpleNamespace(status={}, rank={}, tiles={}, introns=gff,
                            splice_sites=ss, junc_bed=None, params={})
    cfg = _cfg(tmp_path)
    ctrl = Controller(cfg, [], logan=logan)
    assert ctrl._splice_db_path.read_text() == "chr1\t8\t50\t+\n"
    assert ctrl.stranded_introns.counts == {("chr1", 10, 50, "+"): 4}
    assert ctrl._strander.resolve("chr1", 10, 50) == "+"     # cache hit, no genome
    assert ctrl.cumulative_introns.counts == {}                # not merged by default
    with patch("varus.controller.write_hisat2_splice_sites", return_value=1) as m:
        ctrl._rebuild_intron_db()
    m.assert_not_called()                                      # DB already seeded


# ---------------------------------------------------------------------------
# Logan bootstrap, yield scaling, unprocessed discount
# ---------------------------------------------------------------------------

def _logan_ns(**kw):
    base = dict(status={}, rank={}, tiles={}, introns=None, splice_sites=None,
                junc_bed=None, params={}, yield_pct={}, counts={}, acceptance_rate=None)
    base.update(kw)
    return SimpleNamespace(**base)


def test_apply_logan_prior_sets_status_rank_yield(tmp_path: Path):
    cfg = _cfg(tmp_path)
    rng = random.Random(0)
    runs = [RunState.from_record(_rec(a), cfg.batch_size, rng) for a in ("A", "B", "U")]
    logan = _logan_ns(status={"A": "accepted", "B": "accepted"}, rank={"A": 2, "B": 1},
                      tiles={"A": {("c", 0): 1.0}, "B": {("c", 1): 1.0}},
                      yield_pct={"A": 35.0, "B": 250.0})
    kept = apply_logan_prior(runs, logan, batch_size=50_000)
    a, b, u = kept
    assert (a.logan_status, a.logan_rank, a.logan_yield) == ("accepted", 2, 0.35)
    assert b.logan_yield == 1.0                       # clipped
    assert (u.logan_status, u.logan_rank, u.logan_yield) == ("unprocessed", 0, None)


def test_effective_reads_uses_yield_and_unprocessed_weight(tmp_path: Path):
    cfg = _cfg(tmp_path)
    rng = random.Random(0)
    good = RunState.from_record(_rec("G"), cfg.batch_size, rng)
    mixed = RunState.from_record(_rec("M"), cfg.batch_size, rng)
    unproc = RunState.from_record(_rec("U"), cfg.batch_size, rng)
    plain = RunState.from_record(_rec("P"), cfg.batch_size, rng)
    good.logan_status = mixed.logan_status = "accepted"
    good.logan_yield, mixed.logan_yield = 0.8, 0.05
    unproc.logan_status = "unprocessed"
    logan = _logan_ns(counts={"accepted": 3, "rejected": 1, "too_few_contigs": 0},
                      acceptance_rate=0.75)
    ctrl = Controller(cfg, [good, mixed, unproc, plain], logan=logan)
    base = ctrl._effective_reads(plain)
    assert abs(base - 2.0 * cfg.batch_size) < 1e-9          # 100% + 100% priors
    assert abs(ctrl._effective_reads(good) - 0.8 * base) < 1e-9
    assert abs(ctrl._effective_reads(mixed) - 0.05 * base) < 1e-9
    assert abs(ctrl._effective_reads(unproc) - 0.75 * base) < 1e-9
    # explicit weight overrides the acceptance rate; downloaded runs are untouched
    cfg2 = _cfg(tmp_path / "b", logan_unprocessed_weight=1.0)
    ctrl2 = Controller(cfg2, [unproc], logan=logan)
    assert abs(ctrl2._effective_reads(unproc) - base) < 1e-9
    good.times_downloaded = 1
    good.avg_umr_pct, good.avg_spliced_pct = 50.0, 10.0
    assert abs(ctrl._effective_reads(good) - 0.6 * cfg.batch_size) < 1e-9


def test_logan_bootstrap_orders_first_picks(tmp_path: Path, monkeypatch):
    cfg = _cfg(tmp_path, max_batches=5, parallel_downloads=2)
    rng = random.Random(0)
    runs = [RunState.from_record(_rec(a), cfg.batch_size, rng) for a in ("U1", "R2", "R1", "R3")]
    tiles = {a: {("chr1", i): 100} for i, a in enumerate(("U1", "R2", "R1", "R3"))}
    calls = _mock_stack(monkeypatch, tiles)
    logan = _logan_ns(status={"R1": "accepted", "R2": "accepted", "R3": "accepted"},
                      rank={"R1": 1, "R2": 2, "R3": 3},
                      tiles={a: {("chr1", 0): 1.0} for a in ("R1", "R2", "R3")},
                      yield_pct={"R1": 90.0, "R2": 90.0, "R3": 90.0},
                      counts={"accepted": 3, "rejected": 5}, acceptance_rate=3 / 8)
    runs = apply_logan_prior(runs, logan, batch_size=cfg.batch_size)
    ctrl = Controller(cfg, runs, logan=logan)
    # record the pick order itself: download calls come from worker threads,
    # so their order within the in-flight window is racy
    picks = []
    choose = ctrl._choose_next_run_simulated

    def _recording():
        r = choose()
        if r is not None:
            picks.append(r.record.accession)
        return r

    ctrl._choose_next_run_simulated = _recording
    status = ctrl.run()
    assert status == 0
    assert picks[:3] == ["R1", "R2", "R3"], picks
    # disabled bootstrap: cold-start picks are the plain tie-break, not rank order
    cfg2 = _cfg(tmp_path / "nb", max_batches=3, logan_bootstrap=False)
    runs2 = apply_logan_prior(
        [RunState.from_record(_rec(a), cfg2.batch_size, rng) for a in ("U1", "R2", "R1", "R3")],
        logan, batch_size=cfg2.batch_size)
    calls2 = _mock_stack(monkeypatch, tiles)
    ctrl2 = Controller(cfg2, runs2, logan=logan)
    assert not ctrl2._logan_queue


# ---------------------------------------------------------------------------
# Align-ahead: next batch aligns while the current one is scanned/scored
# ---------------------------------------------------------------------------

def _wrap_align_scan(monkeypatch, on_align=None, on_scan=None):
    import varus.controller as vc
    orig_align, orig_scan = vc.align_batch_hisat2, vc.scan_batch_bam
    threads = []

    def align_(*a, **kw):
        threads.append(threading.current_thread().name)
        res = orig_align(*a, **kw)
        if on_align:
            on_align(len(threads))
        return res

    def scan_(*a, **kw):
        if on_scan:
            on_scan()
        return orig_scan(*a, **kw)

    monkeypatch.setattr("varus.controller.align_batch_hisat2", align_)
    monkeypatch.setattr("varus.controller.scan_batch_bam", scan_)
    return threads


def test_align_ahead_overlaps_scan_with_next_alignment(tmp_path: Path, monkeypatch):
    cfg = _cfg(tmp_path, max_batches=4, parallel_downloads=3)
    rng = random.Random(0)
    runs = [RunState.from_record(_rec(a), cfg.batch_size, rng) for a in ("R1", "R2")]
    calls = _mock_stack(monkeypatch, {"R1": {("chr1", 0): 100}, "R2": {("chr1", 1): 100}})
    second_align = threading.Event()
    overlapped = []

    def on_align(n):
        if n == 2:
            second_align.set()

    def on_scan():
        if not overlapped:
            # without align-ahead the 2nd alignment cannot start before this
            # scan returns, so the wait would time out
            overlapped.append(second_align.wait(timeout=5))

    threads = _wrap_align_scan(monkeypatch, on_align, on_scan)
    ctrl = Controller(cfg, runs)
    assert ctrl.run() == 0
    assert overlapped == [True]
    assert all(t.startswith("varus-align") for t in threads), threads
    # the batch aligning ahead counts against max_batches: no overshoot
    assert len(calls["align"]) == 4
    rows = (cfg.outdir / "BatchTimings.tsv").read_text().splitlines()
    assert "t_wait" in rows[0].split("\t") and len(rows) == 5
    assert ctrl._ahead is None and not ctrl._inflight


def test_align_ahead_off_aligns_in_main_thread(tmp_path: Path, monkeypatch):
    cfg = _cfg(tmp_path, max_batches=3, parallel_downloads=3, align_ahead=False)
    rng = random.Random(0)
    runs = [RunState.from_record(_rec(a), cfg.batch_size, rng) for a in ("R1", "R2")]
    _mock_stack(monkeypatch, {"R1": {("chr1", 0): 100}, "R2": {("chr1", 1): 100}})
    threads = _wrap_align_scan(monkeypatch)
    ctrl = Controller(cfg, runs)
    assert ctrl.run() == 0
    assert threads == ["MainThread"] * 3


def test_align_ahead_rejected_batch_is_cleaned_up(tmp_path: Path, monkeypatch):
    cfg = _cfg(tmp_path, max_batches=4, parallel_downloads=3)
    rng = random.Random(0)
    runs = [RunState.from_record(_rec(a), cfg.batch_size, rng) for a in ("BAD", "R2")]
    _mock_stack(monkeypatch, {"BAD": {("chr1", 0): 100}, "R2": {("chr1", 1): 100}},
                uniq={"BAD": 1.0, "R2": 50.0})
    ctrl = Controller(cfg, runs)
    assert ctrl.run() == 0
    bad = next(r for r in runs if r.record.accession == "BAD")
    assert bad.bad_quality and bad.times_downloaded == 0
    assert not (cfg.outdir / "batches" / "BAD").exists() or not any(
        (cfg.outdir / "batches" / "BAD").rglob("Aligned.out.bam"))
    assert ctrl._sim_extra == {}


def test_splice_db_written_atomically(tmp_path: Path):
    cfg = _cfg(tmp_path)
    ctrl = Controller(cfg, [])
    ctrl.stranded_introns = IntronCounts({("chr1", 10, 50, "+"): 1})
    ctrl.cumulative_introns = IntronCounts({("chr1", 10, 50, "."): 1})
    ctrl._db_new_keys = 1
    ctrl._rebuild_intron_db()
    db = cfg.outdir / "intronDB.splice_sites"
    assert db.is_file() and db.read_text().strip()
    assert not (cfg.outdir / "intronDB.splice_sites.tmp").exists()


# ---------------------------------------------------------------------------
# Thread budget: concurrent stages share --threads
# ---------------------------------------------------------------------------

def test_reserve_threads_caps_reservation():
    from varus.align import reserve_threads
    assert reserve_threads(48, 2) == 46
    assert reserve_threads(48, 6) == 42
    assert reserve_threads(48, 30) == 36      # at most a quarter reserved
    assert reserve_threads(4, 3) == 3
    assert reserve_threads(2, 2) == 2         # too few threads to reserve any
    assert reserve_threads(1, 5) == 1


def test_aligner_threads_net_of_concurrent_work(tmp_path: Path):
    from concurrent.futures import Future
    cfg = _cfg(tmp_path, threads=48, parallel_downloads=6)
    ctrl = Controller(cfg, [])
    ctrl._serial = True
    assert ctrl._aligner_threads() == 48                  # nothing runs beside it
    ctrl._serial = False
    assert ctrl._aligner_threads() == 47                  # fastq-dump processes
    ctrl._align_ex = object()
    assert ctrl._aligner_threads() == 46                  # + main thread (align-ahead)
    running = Future()
    ctrl._merge_futures = [running]
    assert ctrl._merge_threads == 4
    assert ctrl._aligner_threads() == 42                  # + rolling merge
    running.set_result(None)
    assert ctrl._aligner_threads() == 46
    ctrl._align_ex = None


def test_align_ahead_passes_budgeted_threads(tmp_path: Path, monkeypatch):
    cfg = _cfg(tmp_path, max_batches=3, parallel_downloads=3, threads=16)
    rng = random.Random(0)
    runs = [RunState.from_record(_rec(a), cfg.batch_size, rng) for a in ("R1", "R2")]
    calls = _mock_stack(monkeypatch, {"R1": {("chr1", 0): 100}, "R2": {("chr1", 1): 100}})
    ctrl = Controller(cfg, runs)
    assert ctrl.run() == 0
    # 16 - main thread - downloads; sort gets a few threads, not 15
    assert {kw["sort_threads"] for _, _, kw in calls["align"]} == {4}


def test_hisat2_sort_threads_flag(tmp_path: Path, monkeypatch):
    seen = []

    class P:
        def __init__(self, cmd, **kw):
            seen.append(cmd)
            self.stdout = type("S", (), {"close": lambda self: None})()
            self.returncode = 0

        def wait(self):
            return 0

        def communicate(self, *a, **k):
            return b"", b""

    monkeypatch.setattr(align.subprocess, "Popen", P)
    monkeypatch.setattr(align, "_require", lambda t: t)
    try:
        align.align_batch_hisat2(tmp_path / "r1.fa", None, index_prefix=tmp_path / "idx",
                                 batch_dir=tmp_path / "b", threads=46, sort_threads=4)
    except Exception:
        pass
    sort = next(c for c in seen if "sort" in c)
    hisat = next(c for c in seen if c[0] == "hisat2")
    assert sort[sort.index("-@") + 1] == "4"
    assert hisat[hisat.index("-p") + 1] == "46"


# ---------------------------------------------------------------------------
# Merged batches: repeated picks of a proven run in one contiguous download
# ---------------------------------------------------------------------------

def _dl_sizes(calls):
    return [(a, (x - n + 1) // 50_000) for a, n, x, _ in calls["download"]]


def test_merge_batches_single_run_merges_contiguous_ranges(tmp_path: Path, monkeypatch):
    cfg = _cfg(tmp_path, max_batches=12, parallel_downloads=2, merge_batches=5)
    rng = random.Random(0)
    runs = [RunState.from_record(_rec("R1", spots=5_000_000), cfg.batch_size, rng)]
    calls = _mock_stack(monkeypatch, {"R1": {("chr1", 0): 100, ("chr1", 1): 50}})
    gate_sizes = []
    import varus.controller as vc
    orig_parse = vc.parse_hisat2_log

    def parse(log_path, batch_size):
        gate_sizes.append(batch_size)
        return orig_parse(log_path, batch_size)

    monkeypatch.setattr("varus.controller.parse_hisat2_log", parse)
    ctrl = Controller(cfg, runs)
    assert ctrl.run() == 0
    sizes = [k for _, k in _dl_sizes(calls)]
    # the first batch is never merged (the run must pass the gate first);
    # downloads run in worker threads, so only the first call's order is fixed
    assert sizes[0] == 1
    assert any(k > 1 for k in sizes) and max(sizes) <= 5
    # the budget is batches, not downloads: exactly max_batches, no overshoot
    assert sum(sizes) == 12 and ctrl.batch_count == 12
    assert runs[0].times_downloaded == 12
    # every download is one contiguous range, and ranges never overlap
    spans = sorted((n, x) for _, n, x, _ in calls["download"])
    assert all(a[1] < b[0] for a, b in zip(spans, spans[1:]))
    # the quality gate divides by the spots actually fetched
    assert sorted(gate_sizes) == sorted(k * 50_000 for k in sizes)
    rows = (cfg.outdir / "BatchTimings.tsv").read_text().splitlines()
    h = rows[0].split("\t")
    assert sum(int(r.split("\t")[h.index("n_batches")]) for r in rows[1:]) == 12


def test_merge_batches_1_and_serial_downloads_never_merge(tmp_path: Path, monkeypatch):
    for kw in (dict(parallel_downloads=2, merge_batches=1),
               dict(parallel_downloads=1, merge_batches=10)):
        cfg = _cfg(tmp_path / str(kw["parallel_downloads"]), max_batches=6, **kw)
        rng = random.Random(0)
        runs = [RunState.from_record(_rec("R1", spots=5_000_000), cfg.batch_size, rng)]
        calls = _mock_stack(monkeypatch, {"R1": {("chr1", 0): 100}})
        assert Controller(cfg, runs).run() == 0
        assert [k for _, k in _dl_sizes(calls)] == [1] * 6, kw


def _merge_ctrl(tmp_path, merge_batches=5):
    cfg = _cfg(tmp_path, merge_batches=merge_batches, max_batches=100)
    rng = random.Random(0)
    good = RunState.from_record(_rec("GOOD", spots=5_000_000), cfg.batch_size, rng)
    other = RunState.from_record(_rec("OTHER", spots=5_000_000), cfg.batch_size, rng)
    good.sigma = list(range(100))           # contiguous order for the test
    ctrl = Controller(cfg, [good, other])
    tiles = [("chr1", i) for i in range(4)]
    ctrl.total_obs = {t: 10 for t in tiles}
    good.times_downloaded = 1
    good.observations = {t: 10 for t in tiles}
    ctrl._estimate_p()                      # builds the tile index and _x_arr
    good.sigma_idx = 1                      # index 0 is the pick being merged
    ctrl._serial = False                    # merging needs pipelined downloads
    return ctrl, good, other


def test_merge_stops_when_another_run_would_win(tmp_path: Path):
    ctrl, good, other = _merge_ctrl(tmp_path)
    exp = {("chr1", 0): 50.0}
    # the other run's (fresh) profit beats GOOD once GOOD's claimed batch counts
    monkeypatch_profit = {id(good): 1.0, id(other): 2.0}
    ctrl._profit = lambda r, extra=None, x_extra=None: monkeypatch_profit[id(r)]
    other.expected_profit = 2.0
    assert ctrl._merge_extra_batches(good, 0, exp) == 1
    assert good.sigma[:3] == [0, 1, 2]      # nothing claimed


def test_merge_claims_only_contiguous_unused_indices(tmp_path: Path):
    ctrl, good, other = _merge_ctrl(tmp_path)
    ctrl._profit = lambda r, extra=None, x_extra=None: 5.0 if r is good else 1.0
    other.expected_profit = 1.0
    good.sigma = [0, 1, 2, 7, 3, 4, 5]      # 3 comes after an unrelated index
    good.sigma_idx = 1
    assert ctrl._merge_extra_batches(good, 0, {("chr1", 0): 1.0}) == 5
    assert good.sigma == [0, 7, 5]          # 1-4 claimed, order of the rest kept
    good.sigma = [0, 2, 3]                  # 1 already used -> no merge
    good.sigma_idx = 1
    assert ctrl._merge_extra_batches(good, 0, {("chr1", 0): 1.0}) == 1


def test_merge_respects_batch_budget_and_needs_a_proven_run(tmp_path: Path):
    ctrl, good, other = _merge_ctrl(tmp_path)
    ctrl._profit = lambda r, extra=None, x_extra=None: 5.0 if r is good else 1.0
    ctrl.batch_count = 97                   # this pick + 2 more fit into 100
    assert ctrl._merge_extra_batches(good, 0, {("chr1", 0): 1.0}) == 3
    ctrl.batch_count = 0
    good.sigma, good.sigma_idx = list(range(100)), 1
    good.times_downloaded = 0               # not yet through the quality gate
    assert ctrl._merge_extra_batches(good, 0, {("chr1", 0): 1.0}) == 1


def test_merged_batches_scanned_in_parallel_singles_in_main_thread(tmp_path: Path, monkeypatch):
    cfg = _cfg(tmp_path, max_batches=12, parallel_downloads=2, merge_batches=5, scan_workers=2)
    rng = random.Random(0)
    runs = [RunState.from_record(_rec("R1", spots=5_000_000), cfg.batch_size, rng)]
    calls = _mock_stack(monkeypatch, {"R1": {("chr1", 0): 100}})
    import varus.controller as vc
    serial, parallel = [], []
    orig_serial, orig_par = vc.scan_batch_bam, vc.scan_batch_bam_parallel

    def s_(bam, tile_size):
        serial.append(bam)
        return orig_serial(bam, tile_size)

    def p_(bam, tile_size, ex, n_parts):
        assert ex is not None and n_parts == 8
        parallel.append(bam)
        return orig_par(bam, tile_size, ex, n_parts)

    monkeypatch.setattr("varus.controller.scan_batch_bam", s_)
    monkeypatch.setattr("varus.controller.scan_batch_bam_parallel", p_)
    ctrl = Controller(cfg, runs)
    assert ctrl.run() == 0
    sizes = [k for _, k in _dl_sizes(calls)]
    assert len(parallel) == sum(1 for k in sizes if k > 1) > 0
    assert len(serial) == sum(1 for k in sizes if k == 1)
    assert ctrl._scan_ex is None                     # pool shut down
    # budget: 2 scan workers reserved instead of the main thread's 1
    ctrl._align_ex, ctrl._scan_ex, ctrl._serial = object(), object(), False
    ctrl.config.threads = 48
    assert ctrl._aligner_threads() == 48 - 2 - 1     # scan workers + downloads
    ctrl._align_ex = ctrl._scan_ex = None


# ---------------------------------------------------------------------------
# Large runlists: lazy sigma and the fresh pool
# ---------------------------------------------------------------------------

def test_lazy_sigma_built_on_first_use_and_deterministic():
    r = RunState.lazy(_rec("L1", spots=260_000), 50_000, seed=42)
    assert r.sigma is None and r.max_batches == 6 and not r.is_exhausted
    n, x = r.next_batch_range(50_000)
    assert r.sigma is not None and sorted(r.sigma) == list(range(6))
    assert r.sigma[-1] == 5                         # short last batch stays last
    assert (n, x) == (r.sigma[0] * 50_000, min(r.sigma[0] * 50_000 + 49_999, 259_999))
    assert RunState.lazy(_rec("L1", spots=260_000), 50_000, seed=42).batch_order() == r.sigma
    r.sigma_idx = 6
    assert r.is_exhausted


def _pool_setup(tmp_path, pool_min, n=40, **kw):
    kw = {"max_batches": 25, "seed": 11, **kw}
    cfg = _cfg(tmp_path / f"p{pool_min}", fresh_pool_min=pool_min, **kw)
    rng = random.Random(3)
    runs = [RunState.from_record(
                _rec(f"SRR{i:03d}", spots=150_000 + 50_000 * (i % 5),
                     avg_len=float(50 + (i * 37) % 200)), cfg.batch_size, rng)
            for i in range(n)]
    tiles = {f"SRR{i:03d}": {("chr1", (i * 3 + k) % 60): 5 + (i + k) % 11 for k in range(6)}
             for i in range(n)}
    return cfg, runs, tiles


def test_fresh_pool_reproduces_serial_pick_sequence(tmp_path: Path, monkeypatch):
    seqs = {}
    for pool_min in (0, 1):
        cfg, runs, tiles = _pool_setup(tmp_path, pool_min)
        calls = _mock_stack(monkeypatch, tiles)
        ctrl = Controller(cfg, runs)
        assert (ctrl._pool_rep is not None) == bool(pool_min)
        assert ctrl.run() == 0
        seqs[pool_min] = [(a, n) for a, n, _, _ in calls["download"]]
        stats = (tmp_path / f"p{pool_min}" / "out" / "RunStatistics.csv").read_text()
        assert len(stats.splitlines()) == 41        # finalize writes every run
    assert len(seqs[0]) == 25
    assert seqs[1] == seqs[0]


def test_fresh_pool_simulated_pick_matches(tmp_path: Path, monkeypatch):
    """Lazy-greedy picks with in-flight extras: pool on == pool off."""
    ctrls = {}
    for pool_min in (0, 1):
        cfg, runs, tiles = _pool_setup(tmp_path, pool_min, max_batches=6)
        _mock_stack(monkeypatch, tiles)
        ctrl = Controller(cfg, runs)
        ctrl.run()                                   # 6 serial batches of state
        ctrl._estimate_p()
        ctrl._calculate_profit()
        ctrls[pool_min] = ctrl
    for extra in ({("chr1", 1): 3000.0, ("chr1", 7): 50.0}, {("chr1", 30): 1e5}, {}):
        picks = []
        for pool_min, ctrl in ctrls.items():
            ctrl._sim_extra = dict(extra)
            got = [ctrl._choose_next_run_simulated().record.accession for _ in range(5)]
            picks.append((got, ctrl.max_profit))
        assert picks[0][0] == picks[1][0]
        assert abs(picks[0][1] - picks[1][1]) < 1e-12


def test_fresh_pool_take_moves_run_and_keeps_order(tmp_path: Path):
    cfg, runs, _ = _pool_setup(tmp_path, 1, n=10)
    runs[4].times_downloaded = 1                     # tracked individually
    ctrl = Controller(cfg, runs)
    assert ctrl._pool_n == 9 and [r.pos for r in ctrl.downloadable] == [4]
    ctrl._estimate_p()
    ctrl._calculate_profit()
    rep = ctrl._pool_rep
    assert rep is runs[0]
    ctrl._pool_take(runs[7])
    ctrl._pool_take(rep)
    assert [r.pos for r in ctrl.downloadable] == [0, 4, 7]
    assert ctrl._pool_rep is runs[1] and ctrl._pool_rep.p_arr is rep.p_arr
    assert ctrl._n_downloadable() == 10 and ctrl._pool_n == 7


def test_fresh_pool_parallel_refill_matches(tmp_path: Path, monkeypatch):
    """_refill with K in-flight picks (pool takes, sim extras): same picks."""
    from concurrent.futures import Future

    class _NoRun:
        def submit(self, fn, *a):
            return Future()

    seqs = []
    for pool_min in (0, 1):
        cfg, runs, tiles = _pool_setup(tmp_path, pool_min, max_batches=8,
                                       merge_batches=4)
        _mock_stack(monkeypatch, tiles)
        ctrl = Controller(cfg, runs)
        ctrl.run()                                   # 8 serial batches of state
        ctrl.config.max_batches = 0
        ctrl._update_downloadable()
        ctrl._estimate_p()
        ctrl._calculate_profit()
        ctrl._serial, ctrl._dl_ex, ctrl._max_inflight = False, _NoRun(), 12
        ctrl._refill()
        seqs.append([(t.run.record.accession, t.n, t.x, t.n_batches) for t in ctrl._inflight])
        assert ctrl._sim_extra
    assert len(seqs[0]) == 12
    assert seqs[1] == seqs[0]


def test_auto_scan_workers_scale_with_threads(tmp_path: Path):
    from varus.controller import auto_scan_workers
    assert [auto_scan_workers(t) for t in (1, 4, 8, 15, 16, 24, 32, 48, 256)] == \
        [0, 0, 0, 0, 2, 3, 4, 4, 4]
    cfg = _cfg(tmp_path, threads=8, merge_batches=5, parallel_downloads=2)
    assert cfg.scan_workers is None
    Controller(cfg, [])
    assert cfg.scan_workers == 0

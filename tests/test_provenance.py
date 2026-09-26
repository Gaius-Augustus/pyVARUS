"""Tests for varus.provenance and varus.replay: manifest, splice-DB log, replay.

Download, alignment and merge are mocked; the test checks that a replay
fetches exactly the batches of the original run and aligns each against the
same splice-site DB content.
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

import varus.controller as vc
from varus import download
from varus.controller import Controller, RunState
from varus.download import BatchPaths
from varus.introns import IntronCounts
from varus.provenance import (
    MANIFEST_NAME, SPLICE_LOG_NAME, SpliceDBLog, SpliceDBSnapshots,
    read_manifest, write_manifest,
)
from varus.replay import EXIT_INCOMPLETE, ReplayConfig, replay
from varus.tiles import BAMStats

from tests.conftest import requires_pysam
from tests.test_speedups import _cfg, _mock_stack, _rec, _write_genome


# ---------------------------------------------------------------------------
# Splice-DB log
# ---------------------------------------------------------------------------

def test_splice_log_roundtrip(tmp_path: Path):
    db = tmp_path / "db"
    lg = SpliceDBLog(tmp_path / "log.gz", "hisat2")
    assert lg.version == 0
    db.write_text("chr1\t10\t50\t+\n")
    assert lg.record(db) == 1
    assert lg.record(db) == 1                       # unchanged -> same version
    db.write_text("chr1\t10\t50\t+\nchr1\t9\t40\t-\nchr2\t5\t9\t+\n")
    assert lg.record(db) == 2
    db.write_text("chr1\t9\t40\t-\nchr2\t5\t9\t+\n")  # a line removed (BED score change)
    assert lg.record(db) == 3

    snaps = SpliceDBSnapshots(tmp_path / "log.gz")
    assert snaps.fmt == "hisat2" and snaps.max_version == 3
    assert snaps.lines(1) == ["chr1\t10\t50\t+"]
    # numeric order, like the writers
    assert snaps.lines(2) == ["chr1\t9\t40\t-", "chr1\t10\t50\t+", "chr2\t5\t9\t+"]
    assert snaps.lines(3) == ["chr1\t9\t40\t-", "chr2\t5\t9\t+"]
    assert snaps.lines(1) == ["chr1\t10\t50\t+"]     # rewinds
    assert snaps.write(0, tmp_path / "none") is None
    with pytest.raises(ValueError):
        snaps.lines(4)


def test_manifest_roundtrip(tmp_path: Path):
    rows = [{"batch": 1, "accession": "SRR1", "n": 0, "x": 49_999, "paired": 1,
             "platform": "", "preset": "", "n_batches": 1, "db_version": 0,
             "align_threads": 4, "uniq_pct": "81.2500"}]
    p = write_manifest(tmp_path / MANIFEST_NAME, {"genome_md5": "abc", "command": "a\tb"},
                       rows)
    header, got = read_manifest(p)
    assert header["genome_md5"] == "abc" and header["command"] == "a b"
    assert got == [{"batch": 1, "accession": "SRR1", "n": 0, "x": 49_999,
                    "paired": True, "platform": "", "preset": "", "n_batches": 1,
                    "db_version": 0, "align_threads": 4, "uniq_pct": 81.25}]


# ---------------------------------------------------------------------------
# Controller writes what replay needs
# ---------------------------------------------------------------------------

def _mock_growing_db(monkeypatch, n_runs=3):
    """Mock stack whose every batch reports one new intron, so the DB grows."""
    calls = _mock_stack(monkeypatch, {f"R{i}": {("chr1", i): 100} for i in range(n_runs)})
    seen = []

    def fake_scan(bam, tile_size):
        acc = bam.parent.parent.name
        k = len(seen)
        seen.append(acc)
        stats = BAMStats(umr_counts={("chr1", int(acc[1:])): 100}, n_reads=100,
                         n_spliced=0)
        # _write_genome's canonical chr1 intron, plus a new junction each batch
        return stats, IntronCounts({("chr1", 11, 50, "."): 1,
                                    ("chr1", 11, 60 + k, "."): 1})

    monkeypatch.setattr("varus.controller.scan_batch_bam", fake_scan)
    monkeypatch.setattr("varus.controller.scan_batch_bam_parallel",
                        lambda bam, tile_size, ex, n_parts: fake_scan(bam, tile_size))

    db_seen = {}
    orig_align = vc.align_batch_hisat2

    def align(r1, r2, *, index_prefix, batch_dir, threads, intron_db=None, **kw):
        db_seen[batch_dir.name + batch_dir.parent.name] = (
            intron_db.read_text() if intron_db is not None else None, threads)
        return orig_align(r1, r2, index_prefix=index_prefix, batch_dir=batch_dir,
                          threads=threads, intron_db=intron_db, **kw)

    monkeypatch.setattr("varus.controller.align_batch_hisat2", align)
    return calls, db_seen


@requires_pysam
@pytest.mark.parametrize("kw", [dict(parallel_downloads=1),
                                dict(parallel_downloads=3, merge_batches=3),
                                dict(parallel_downloads=1, bootstrap_all=True)])
def test_run_writes_manifest_and_db_log(tmp_path: Path, monkeypatch, kw):
    cfg = _cfg(tmp_path, max_batches=8, **kw)
    _write_genome(cfg.genome)
    rng = random.Random(0)
    runs = [RunState.from_record(_rec(f"R{i}", spots=5_000_000), cfg.batch_size, rng)
            for i in range(3)]
    calls, db_seen = _mock_growing_db(monkeypatch)
    ctrl = Controller(cfg, runs)
    assert ctrl.run() == 0

    header, rows = read_manifest(cfg.outdir / MANIFEST_NAME)
    assert sum(r["n_batches"] for r in rows) == ctrl.batch_count == int(header["batches"])
    assert len(rows) == len(calls["download"])
    assert sorted((r["accession"], r["n"], r["x"]) for r in rows) == sorted(calls["download"])
    assert header["genome_md5"] and header["aligner"] == "hisat2"
    # the first batch has no DB; later ones see a growing one
    assert rows[0]["db_version"] == 0
    assert rows[-1]["db_version"] > 0
    versions = [r["db_version"] for r in rows]
    assert versions == sorted(versions)
    # each batch's DB content is exactly the logged version
    snaps = SpliceDBSnapshots(cfg.outdir / SPLICE_LOG_NAME)
    for r in rows:
        key = f"N{r['n']}X{r['x']}{r['accession']}"
        want = None if r["db_version"] == 0 else "".join(
            ln + "\n" for ln in snaps.lines(r["db_version"]))
        assert db_seen[key] == (want, r["align_threads"]), r
    # the pinned DB links are cleaned up with the batch
    assert not list((cfg.outdir).glob("batches/**/intronDB.splice_sites"))


# ---------------------------------------------------------------------------
# Replay
# ---------------------------------------------------------------------------

def _mock_replay(monkeypatch, fail=()):
    got = {"download": [], "db": {}, "merged": []}

    def fake_download(accession, n, x, paired, outdir, **kw):
        if accession in fail:
            raise RuntimeError("gone")
        bdir = download.batch_dir_for(outdir, accession, n, x)
        bdir.mkdir(parents=True, exist_ok=True)
        r1 = bdir / f"{accession}.fasta"
        r1.write_text(">r\nACGT\n")
        got["download"].append((accession, n, x))
        return BatchPaths(r1=r1, r2=None, batch_dir=bdir)

    def fake_align(r1, r2, *, index_prefix, batch_dir, threads, intron_db=None, **kw):
        got["db"][f"{batch_dir.name}{batch_dir.parent.name}"] = (
            intron_db.read_text() if intron_db is not None else None, threads)
        bam = batch_dir / "Aligned.out.bam"
        bam.write_bytes(b"BAM")
        logp = batch_dir / "Log.final.out"
        logp.write_text("x")
        from varus.align import AlignmentResult
        return AlignmentResult(bam=bam, log=logp)

    def fake_merge(bams, out, threads=4, compression=None, **kw):
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(b"MERGED")
        got["merged"].append((list(bams), out))
        return out

    monkeypatch.setattr("varus.replay.download_batch", fake_download)
    monkeypatch.setattr("varus.replay.align_batch_hisat2", fake_align)
    monkeypatch.setattr("varus.replay.parse_hisat2_log",
                        lambda log_path, batch_size: {"uniq_pct": 50.0})
    monkeypatch.setattr("varus.replay.merge_bams", fake_merge)
    return got


@requires_pysam
def test_replay_fetches_same_batches_with_same_db(tmp_path: Path, monkeypatch):
    cfg = _cfg(tmp_path / "orig", max_batches=8, parallel_downloads=3, merge_batches=3)
    cfg.outdir.parent.mkdir(parents=True, exist_ok=True)
    _write_genome(cfg.genome)
    rng = random.Random(0)
    runs = [RunState.from_record(_rec(f"R{i}", spots=5_000_000), cfg.batch_size, rng)
            for i in range(3)]
    calls, db_seen = _mock_growing_db(monkeypatch)
    assert Controller(cfg, runs).run() == 0
    manifest = cfg.outdir / MANIFEST_NAME

    got = _mock_replay(monkeypatch)
    out = tmp_path / "replayed"
    rc = replay(ReplayConfig(manifest=manifest, genome=cfg.genome, index=cfg.index_prefix,
                             outdir=out, parallel_downloads=2))
    assert rc == 0
    assert sorted(got["download"]) == sorted(calls["download"])
    assert got["db"] == db_seen
    assert (out / "VARUS.bam").read_bytes() == b"MERGED"
    assert not (out / "replay").exists()          # batch files cleaned up

    # a run that SRA no longer serves: incomplete BAM, list of what is missing
    _, rows = read_manifest(manifest)
    gone = rows[-1]["accession"]
    _mock_replay(monkeypatch, fail={gone})
    out2 = tmp_path / "replayed2"
    rc = replay(ReplayConfig(manifest=manifest, genome=cfg.genome, index=cfg.index_prefix,
                             outdir=out2))
    assert rc == EXIT_INCOMPLETE
    assert not (out2 / "VARUS.bam").exists()
    assert (out2 / "VARUS.incomplete.bam").is_file()
    assert gone in (out2 / "replay_missing.tsv").read_text()


def test_replay_refuses_other_genome(tmp_path: Path, monkeypatch):
    g = tmp_path / "g.fa"
    g.write_text(">chr1\nACGT\n")
    manifest = write_manifest(
        tmp_path / MANIFEST_NAME, {"genome": "g.fa", "genome_md5": "0" * 32},
        [{"batch": 1, "accession": "SRR1", "n": 0, "x": 9, "paired": 0, "platform": "",
          "preset": "", "n_batches": 1, "db_version": 0, "align_threads": 4, "uniq_pct": "50.0000"}])
    _mock_replay(monkeypatch)
    with pytest.raises(SystemExit, match="MD5"):
        replay(ReplayConfig(manifest=manifest, genome=g, index=tmp_path / "idx",
                            outdir=tmp_path / "o"))
    rc = replay(ReplayConfig(manifest=manifest, genome=g, index=tmp_path / "idx",
                             outdir=tmp_path / "o", skip_genome_check=True))
    assert rc == 0


def test_replay_needs_db_log(tmp_path: Path):
    manifest = write_manifest(
        tmp_path / MANIFEST_NAME, {},
        [{"batch": 1, "accession": "SRR1", "n": 0, "x": 9, "paired": 0, "platform": "",
          "preset": "", "n_batches": 1, "db_version": 2, "align_threads": 4, "uniq_pct": "50.0000"}])
    with pytest.raises(SystemExit, match="splice-db-log"):
        replay(ReplayConfig(manifest=manifest, genome=tmp_path / "g.fa",
                            index=tmp_path / "idx", outdir=tmp_path / "o"))



def test_replay_longreads_uses_minimap2_preset_and_bed(tmp_path: Path, monkeypatch):
    db = tmp_path / "db"
    lg = SpliceDBLog(tmp_path / SPLICE_LOG_NAME, "bed12")
    db.write_text("chr1\t8\t51\t.\t3\t+\n")
    assert lg.record(db) == 1
    manifest = write_manifest(
        tmp_path / MANIFEST_NAME, {"mode": "longreads", "min_mapq": "1"},
        [{"batch": 1, "accession": "SRR9", "n": 0, "x": 1999, "paired": 0,
          "platform": "OXFORD_NANOPORE", "preset": "ont", "n_batches": 1,
          "db_version": 1, "align_threads": 4, "uniq_pct": "70.0000"}])
    got = _mock_replay(monkeypatch)
    seen = {}

    def fake_mm2(reads, *, index, batch_dir, threads, preset, junc_bed=None, **kw):
        seen.update(preset=preset, bed=junc_bed.read_text(), index=index)
        bam = batch_dir / "Aligned.out.bam"
        bam.write_bytes(b"BAM")
        from varus.align import AlignmentResult
        return AlignmentResult(bam=bam, log=batch_dir / "Log.minimap2.err")

    monkeypatch.setattr("varus.replay.align_batch_minimap2", fake_mm2)
    monkeypatch.setattr("varus.replay.count_minimap2_quality",
                        lambda bam, min_mapq: {"uniq_pct": 70.0})
    rc = replay(ReplayConfig(manifest=manifest, genome=tmp_path / "g.fa",
                             index=tmp_path / "mm2idx.mmi", outdir=tmp_path / "o"))
    assert rc == 0
    assert seen == {"preset": "ont", "bed": "chr1\t8\t51\t.\t3\t+\n",
                    "index": tmp_path / "mm2idx.mmi"}
    assert got["download"] == [("SRR9", 0, 1999)]

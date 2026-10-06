"""Tests for sampling from local gzipped FASTQ files (``varus run --fastq``).

Unit tests for the gzip index and record extraction, the wiring into the
controller/CLI/replay with mocked aligners, and one small end-to-end run
with real HISAT2 + samtools when they are on PATH.
"""

from __future__ import annotations

import gzip
import io
import random
import shutil
from pathlib import Path

import pytest

from varus import cli
from varus.controller import Controller, RunState, load_local_runs
from varus.localreads import (
    GzipFastqIndex, LocalReadsError, extract_local_batch, index_local_runs,
    parse_fastq_specs, run_name,
)
from varus.provenance import MANIFEST_NAME, read_manifest

from tests.conftest import requires_pysam
from tests.test_speedups import _cfg, _mock_stack, _write_genome


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _reads(n: int, rng: random.Random, lmin=30, lmax=120, prefix="r"):
    out = []
    for i in range(n):
        L = rng.randint(lmin, lmax)
        out.append((f"{prefix}{i}".encode(),
                    "".join(rng.choice("ACGT") for _ in range(L)).encode()))
    return out


def _fastq_bytes(reads, comment=b" extra words", crlf=False) -> bytes:
    nl = b"\r\n" if crlf else b"\n"
    buf = io.BytesIO()
    for name, seq in reads:
        buf.write(b"@" + name + comment + nl + seq + nl + b"+" + nl + b"I" * len(seq) + nl)
    return buf.getvalue()


def _write_gz(path: Path, raw: bytes, members: int = 1, rng=None) -> Path:
    """gzip ``raw`` into ``path`` as one or several concatenated members."""
    if members <= 1:
        path.write_bytes(gzip.compress(raw))
        return path
    rng = rng or random.Random(0)
    cuts = sorted(rng.sample(range(1, len(raw)), members - 1))
    parts, prev = [], 0
    for c in cuts + [len(raw)]:
        parts.append(gzip.compress(raw[prev:c]))
        prev = c
    parts.append(gzip.compress(b""))        # bgzip-style empty EOF member
    path.write_bytes(b"".join(parts))
    return path


# Small chunks/snapshots/strides so a few-MB test file exercises every path:
# resume from a snapshot, skip records after a stored start, member boundaries.
_SMALL = dict(stride=7, snapshot_bytes=20_000, chunk_bytes=4096)


# ---------------------------------------------------------------------------
# GzipFastqIndex
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("members,crlf", [(1, False), (5, False), (3, True)])
def test_index_counts_and_extracts_exactly(tmp_path: Path, members, crlf):
    rng = random.Random(members)
    reads = _reads(3000, rng)
    fq = _write_gz(tmp_path / "lib.fq.gz", _fastq_bytes(reads, crlf=crlf), members, rng)
    idx = GzipFastqIndex(fq, **_SMALL)
    assert idx.n_records == 3000
    assert idx.n_bases == sum(len(s) for _, s in reads)
    assert len(idx._snapshots) > 2 and len(idx._starts) == (3000 - 1) // 7 + 1
    for _ in range(100):
        n = rng.randrange(3000)
        x = min(2999, n + rng.randrange(0, 300))
        assert list(idx.records(n, x)) == reads[n:x + 1]
    # the last batch is capped at the file end; the very first/last spot work
    assert list(idx.records(2999, 5000)) == reads[-1:]
    assert list(idx.records(0, 0)) == reads[:1]
    with pytest.raises(ValueError):
        list(idx.records(3000, 3001))


def test_index_without_trailing_newline(tmp_path: Path):
    reads = _reads(10, random.Random(1))
    raw = _fastq_bytes(reads, comment=b"")[:-1]
    idx = GzipFastqIndex(_write_gz(tmp_path / "a.fastq.gz", raw))
    assert idx.n_records == 10 and idx.n_bases == sum(len(s) for _, s in reads)
    assert list(idx.records(8, 9)) == reads[8:]


def test_write_fasta(tmp_path: Path):
    reads = _reads(50, random.Random(2))
    idx = GzipFastqIndex(_write_gz(tmp_path / "a.fq.gz", _fastq_bytes(reads)), **_SMALL)
    out = tmp_path / "b.fasta"
    assert idx.write_fasta(10, 14, out) == 5
    text = out.read_text()
    assert text == "".join(f">{n.decode()}\n{s.decode()}\n" for n, s in reads[10:15])
    assert not out.with_name("b.fasta.tmp").exists()


@pytest.mark.parametrize("make,match", [
    (lambda raw: raw, "not gzip-compressed"),
    (lambda raw: gzip.compress(raw[:-10]), "quality and sequence"),      # truncated qualities
    (lambda raw: gzip.compress(raw + b"@orphan\nACGT\n"), "not a multiple of 4"),
    (lambda raw: gzip.compress(raw.replace(b"@", b">")), "does not start with '@'"),
    (lambda raw: gzip.compress(raw)[:-30], "corrupt"),
])
def test_bad_files_are_rejected(tmp_path: Path, make, match):
    raw = _fastq_bytes(_reads(20, random.Random(3)), comment=b"")
    fq = tmp_path / "bad.fq.gz"
    fq.write_bytes(make(raw))
    with pytest.raises(LocalReadsError, match=match):
        idx = GzipFastqIndex(fq)
        list(idx.records(0, 19))


def test_wrong_suffix_is_rejected(tmp_path: Path):
    raw = _fastq_bytes(_reads(5, random.Random(4)))
    for name in ("reads.fastq", "reads.fq", "reads.fasta.gz", "reads.gz"):
        p = tmp_path / name
        p.write_bytes(gzip.compress(raw))
        with pytest.raises(LocalReadsError, match="fastq.gz"):
            GzipFastqIndex(p)
    with pytest.raises(LocalReadsError, match="not found"):
        GzipFastqIndex(tmp_path / "missing.fq.gz")


def test_empty_file_is_rejected(tmp_path: Path):
    p = tmp_path / "empty.fq.gz"
    p.write_bytes(gzip.compress(b""))
    with pytest.raises(LocalReadsError, match="no reads"):
        GzipFastqIndex(p)


# ---------------------------------------------------------------------------
# Naming and --fastq parsing
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("r1,r2,want", [
    ("reads.fq.gz", None, "reads"),
    ("lib_1.fq.gz", "lib_2.fq.gz", "lib"),
    ("S1_R1_001.fastq.gz", "S1_R2_001.fastq.gz", "S1"),
    ("liver.R1.fq.gz", "liver.R2.fq.gz", "liver"),
    ("a_fwd.fq.gz", "a_rev.fq.gz", "a"),
    ("x.fq.gz", "y.fq.gz", "x"),
    ("1.fq.gz", "2.fq.gz", "1"),
    ("my sample (2).fq.gz", None, "my_sample_2"),
])
def test_run_name(r1, r2, want):
    assert run_name(Path("/data") / r1, Path("/data") / r2 if r2 else None) == want


def test_parse_fastq_specs():
    specs = parse_fastq_specs(["a.fq.gz", "b_1.fq.gz,b_2.fq.gz"])
    assert specs == [(Path("a.fq.gz"), None), (Path("b_1.fq.gz"), Path("b_2.fq.gz"))]
    with pytest.raises(LocalReadsError, match="one file"):
        parse_fastq_specs(["a,b,c"])
    with pytest.raises(LocalReadsError, match="same file"):
        parse_fastq_specs(["a.fq.gz,a.fq.gz"])


def test_index_local_runs_paired_and_errors(tmp_path: Path):
    rng = random.Random(5)
    r1 = _reads(40, rng, prefix="p")
    r2 = _reads(40, rng, prefix="p")
    _write_gz(tmp_path / "lib_1.fq.gz", _fastq_bytes(r1))
    _write_gz(tmp_path / "lib_2.fq.gz", _fastq_bytes(r2))
    _write_gz(tmp_path / "se.fq.gz", _fastq_bytes(_reads(13, rng)))
    runs = index_local_runs([(tmp_path / "lib_1.fq.gz", tmp_path / "lib_2.fq.gz"),
                             (tmp_path / "se.fq.gz", None)], threads=3, **_SMALL)
    assert [r.name for r in runs] == ["lib", "se"]
    assert runs[0].paired and runs[0].total_spots == 40
    assert runs[0].total_bases == sum(len(s) for _, s in r1 + r2)
    assert runs[0].avg_len == pytest.approx(runs[0].total_bases / 40)
    rec = runs[0].record("ILLUMINA")
    assert rec.accession == "lib" and rec.paired and rec.platform == "ILLUMINA"
    assert not runs[1].paired and runs[1].total_spots == 13

    # mates of different length
    _write_gz(tmp_path / "odd_2.fq.gz", _fastq_bytes(_reads(39, rng)))
    with pytest.raises(LocalReadsError, match="different numbers of reads"):
        index_local_runs([(tmp_path / "lib_1.fq.gz", tmp_path / "odd_2.fq.gz")])
    # two libraries with the same name
    (tmp_path / "d").mkdir()
    shutil.copy(tmp_path / "se.fq.gz", tmp_path / "d" / "se.fq.gz")
    with pytest.raises(LocalReadsError, match="same name"):
        index_local_runs([(tmp_path / "se.fq.gz", None), (tmp_path / "d" / "se.fq.gz", None)])
    # explicit names (replay)
    runs = index_local_runs([(tmp_path / "se.fq.gz", None)], names=["X1"])
    assert runs[0].name == "X1"


def test_extract_local_batch_layout(tmp_path: Path):
    rng = random.Random(6)
    r1, r2 = _reads(30, rng, prefix="a"), _reads(30, rng, prefix="b")
    _write_gz(tmp_path / "lib_1.fq.gz", _fastq_bytes(r1))
    _write_gz(tmp_path / "lib_2.fq.gz", _fastq_bytes(r2))
    _write_gz(tmp_path / "se.fq.gz", _fastq_bytes(r1))
    pe, se = index_local_runs([(tmp_path / "lib_1.fq.gz", tmp_path / "lib_2.fq.gz"),
                               (tmp_path / "se.fq.gz", None)], **_SMALL)
    out = tmp_path / "out"
    paths = extract_local_batch(pe, 10, 19, out)
    assert paths.batch_dir == out / "batches" / "lib" / "N10X19"
    assert paths.r1.name == "lib_1.fasta" and paths.r2.name == "lib_2.fasta"
    assert paths.r1.read_text().count(">") == 10
    assert paths.r2.read_text().splitlines()[1] == r2[10][1].decode()
    paths = extract_local_batch(se, 20, 99, out)      # capped at the end
    assert paths.r2 is None and paths.r1 == out / "batches" / "se" / "N20X99" / "se.fasta"
    assert paths.r1.read_text().count(">") == 10
    with pytest.raises(RuntimeError, match="extracting spots"):
        extract_local_batch(se, 30, 39, out)


# ---------------------------------------------------------------------------
# Controller: a local run is sampled like an SRA run
# ---------------------------------------------------------------------------

@requires_pysam
def test_controller_samples_local_run_batches(tmp_path: Path, monkeypatch):
    rng = random.Random(7)
    reads = _reads(1000, rng)
    _write_gz(tmp_path / "lib.fq.gz", _fastq_bytes(reads), members=3, rng=rng)
    local = index_local_runs([(tmp_path / "lib.fq.gz", None)], **_SMALL)
    cfg = _cfg(tmp_path, batch_size=100, max_batches=4, keep_batches=True)
    _write_genome(cfg.genome)
    runs = load_local_runs(local, cfg.batch_size, random.Random(0))
    assert len(runs) == 1 and runs[0].max_batches == 10 and runs[0].local is local[0]
    calls = _mock_stack(monkeypatch, {"lib": {("chr1", 0): 5}})
    ctrl = Controller(cfg, runs)
    assert ctrl.run() == 0
    assert ctrl.batch_count == 4
    assert calls["download"] == []            # fastq-dump was never called
    header, rows = read_manifest(cfg.outdir / MANIFEST_NAME)
    assert [r["accession"] for r in rows] == ["lib"] * 4
    assert header["local_fastq.lib"] == str((tmp_path / "lib.fq.gz").resolve())
    # every batch holds exactly the spots of its range, in file order
    ns = sorted(r["n"] for r in rows)
    assert len(set(ns)) == 4 and all(n % 100 == 0 for n in ns)
    for r in rows:
        fa = cfg.outdir / "batches" / "lib" / f"N{r['n']}X{r['x']}" / "lib.fasta"
        names = [ln[1:] for ln in fa.read_text().splitlines() if ln.startswith(">")]
        assert names == [nm.decode() for nm, _ in reads[r["n"]:r["x"] + 1]]
    stats = (cfg.outdir / "RunStatistics.csv").read_text()
    assert "lib;4;" in stats


@requires_pysam
def test_controller_mixes_sra_and_local_runs(tmp_path: Path, monkeypatch):
    from tests.test_speedups import _rec
    rng = random.Random(8)
    _write_gz(tmp_path / "lib.fq.gz", _fastq_bytes(_reads(500, rng)))
    local = index_local_runs([(tmp_path / "lib.fq.gz", None)])
    cfg = _cfg(tmp_path, batch_size=100, max_batches=6, parallel_downloads=3, merge_batches=1)
    _write_genome(cfg.genome)
    runs = [RunState.from_record(_rec("SRR1", spots=1000), cfg.batch_size, rng)]
    runs += load_local_runs(local, cfg.batch_size, rng)
    # the SRA run yields new tiles every time, the local run only one
    calls = _mock_stack(monkeypatch, {"SRR1": {("chr1", 0): 3, ("chr2", 0): 3},
                                      "lib": {("chr1", 0): 3}})
    ctrl = Controller(cfg, runs)
    assert ctrl.run() == 0
    _, rows = read_manifest(cfg.outdir / MANIFEST_NAME)
    accs = {r["accession"] for r in rows}
    assert accs == {"SRR1", "lib"}
    assert len(calls["download"]) == sum(1 for r in rows if r["accession"] == "SRR1")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def test_run_parser_accepts_fastq(tmp_path):
    args = cli.build_parser().parse_args(
        ["run", "Foo bar", "g.fa", "--index", "idx", "--fastq", "a.fq.gz",
         "--fastq", "b_1.fq.gz,b_2.fq.gz", "--fastq-platform", "OXFORD_NANOPORE"])
    assert args.runlist is None
    assert args.fastq == ["a.fq.gz", "b_1.fq.gz,b_2.fq.gz"]
    assert args.fastq_platform == "OXFORD_NANOPORE"
    args = cli.build_parser().parse_args(["replay", "m.tsv", "g.fa", "--index", "i",
                                          "--fastq-dir", str(tmp_path)])
    assert args.fastq_dir == tmp_path


def test_run_requires_runlist_or_fastq(tmp_path, monkeypatch):
    (tmp_path / "g.fa").write_text(">c\nACGT\n")
    with pytest.raises(SystemExit, match="--runlist.*--fastq"):
        cli.main(["run", "Foo bar", str(tmp_path / "g.fa"), "--index", "i", "--no-logan"])


def _fake_controller(monkeypatch):
    import varus.controller as ctl
    seen = {}

    class FakeController:
        def __init__(self, cfg, runs, logan=None):
            seen["runs"] = runs
            seen["cfg"] = cfg
            seen["logan"] = logan

        def run(self):
            return 0

    monkeypatch.setattr(ctl, "Controller", FakeController)
    return seen


def test_run_with_fastq_only_skips_logan_and_builds_local_runs(tmp_path, monkeypatch):
    rng = random.Random(9)
    _write_gz(tmp_path / "lib_1.fq.gz", _fastq_bytes(_reads(60, rng)))
    _write_gz(tmp_path / "lib_2.fq.gz", _fastq_bytes(_reads(60, rng)))
    _write_gz(tmp_path / "se.fq.gz", _fastq_bytes(_reads(25, rng)))
    (tmp_path / "g.fa").write_text(">c\nACGT\n")
    monkeypatch.setattr(cli, "logan_prescreen",
                        lambda a, c: pytest.fail("Logan pre-screen must not run"))
    seen = _fake_controller(monkeypatch)
    rc = cli.main(["run", "Foo bar", str(tmp_path / "g.fa"), "--index", "idx",
                   "--outdir", str(tmp_path), "--batch-size", "10", "--seed", "1",
                   "--fastq", f"{tmp_path / 'lib_1.fq.gz'},{tmp_path / 'lib_2.fq.gz'}",
                   "--fastq", str(tmp_path / "se.fq.gz")])
    assert rc == 0
    runs = seen["runs"]
    assert [r.record.accession for r in runs] == ["lib", "se"]
    assert all(r.local is not None for r in runs)
    assert runs[0].record.paired and runs[0].max_batches == 6
    assert not runs[1].record.paired and runs[1].max_batches == 3
    assert seen["logan"] is None


def test_run_with_runlist_and_fastq_appends_local_runs(tmp_path, monkeypatch):
    import varus.controller as ctl
    from varus.runlist import RunRecord
    rng = random.Random(10)
    _write_gz(tmp_path / "se.fq.gz", _fastq_bytes(_reads(25, rng)))
    (tmp_path / "g.fa").write_text(">c\nACGT\n")
    rl = tmp_path / "Runlist.tsv"
    rl.write_text("")

    def fake_load_runs(path, batch_size, rng):
        rec = RunRecord(accession="SRR9", total_spots=100, total_bases=10_000, avg_len=100.0,
                        paired=False, colorspace=False, platform="ILLUMINA")
        return [ctl.RunState.from_record(rec, batch_size, rng)]

    monkeypatch.setattr(ctl, "load_runs", fake_load_runs)
    seen = _fake_controller(monkeypatch)
    rc = cli.main(["run", "Foo bar", str(tmp_path / "g.fa"), "--index", "idx", "--no-logan",
                   "--runlist", str(rl), "--outdir", str(tmp_path), "--batch-size", "10",
                   "--fastq", str(tmp_path / "se.fq.gz"), "--longreads",
                   "--fastq-platform", "OXFORD_NANOPORE"])
    assert rc == 0
    accs = [r.record.accession for r in seen["runs"]]
    assert accs == ["SRR9", "se"]
    assert seen["runs"][1].record.platform == "OXFORD_NANOPORE"


def test_run_reports_bad_fastq(tmp_path, monkeypatch):
    (tmp_path / "g.fa").write_text(">c\nACGT\n")
    (tmp_path / "plain.fq.gz").write_bytes(b"@r\nACGT\n+\nIIII\n")
    _fake_controller(monkeypatch)
    with pytest.raises(SystemExit, match="not gzip-compressed"):
        cli.main(["run", "Foo bar", str(tmp_path / "g.fa"), "--index", "idx", "--no-logan",
                  "--outdir", str(tmp_path), "--fastq", str(tmp_path / "plain.fq.gz")])
    with pytest.raises(SystemExit, match="one file"):
        cli.main(["run", "Foo bar", str(tmp_path / "g.fa"), "--index", "idx", "--no-logan",
                  "--outdir", str(tmp_path), "--fastq", "a.fq.gz,b.fq.gz,c.fq.gz"])


# ---------------------------------------------------------------------------
# Replay
# ---------------------------------------------------------------------------

@requires_pysam
def test_replay_extracts_local_batches(tmp_path: Path, monkeypatch):
    from tests.test_provenance import _mock_replay
    from varus.replay import ReplayConfig, replay
    rng = random.Random(11)
    reads = _reads(300, rng)
    data = tmp_path / "data"
    data.mkdir()
    _write_gz(data / "lib.fq.gz", _fastq_bytes(reads))
    local = index_local_runs([(data / "lib.fq.gz", None)], **_SMALL)
    cfg = _cfg(tmp_path, batch_size=50, max_batches=3, outdir=tmp_path / "run")
    _write_genome(cfg.genome)
    _mock_stack(monkeypatch, {"lib": {("chr1", 0): 5}})
    ctrl = Controller(cfg, load_local_runs(local, cfg.batch_size, random.Random(0)))
    assert ctrl.run() == 0
    header, rows = read_manifest(cfg.outdir / MANIFEST_NAME)
    assert len(rows) == 3

    got = _mock_replay(monkeypatch)
    out = tmp_path / "replay"
    rc = replay(ReplayConfig(manifest=cfg.outdir / MANIFEST_NAME, genome=cfg.genome,
                             index=tmp_path / "idx", outdir=out, keep_batches=True))
    assert rc == 0
    assert got["download"] == []                       # no SRA download
    for r in rows:
        fa = out / "replay" / "batches" / "lib" / f"N{r['n']}X{r['x']}" / "lib.fasta"
        names = [ln[1:] for ln in fa.read_text().splitlines() if ln.startswith(">")]
        assert names == [nm.decode() for nm, _ in reads[r["n"]:r["x"] + 1]]

    # file moved: --fastq-dir finds it by name, without it replay stops
    moved = tmp_path / "moved"
    moved.mkdir()
    shutil.move(str(data / "lib.fq.gz"), str(moved / "lib.fq.gz"))
    with pytest.raises(SystemExit, match="--fastq-dir"):
        replay(ReplayConfig(manifest=cfg.outdir / MANIFEST_NAME, genome=cfg.genome,
                            index=tmp_path / "idx", outdir=tmp_path / "replay2"))
    rc = replay(ReplayConfig(manifest=cfg.outdir / MANIFEST_NAME, genome=cfg.genome,
                             index=tmp_path / "idx", outdir=tmp_path / "replay3",
                             fastq_dir=moved))
    assert rc == 0


# ---------------------------------------------------------------------------
# End to end with the real aligner (small data)
# ---------------------------------------------------------------------------

_TOOLS = all(shutil.which(t) for t in ("hisat2", "hisat2-build", "samtools"))


def _random_genome(rng: random.Random, lengths) -> dict:
    return {f"chr{i + 1}": "".join(rng.choice("ACGT") for _ in range(L))
            for i, L in enumerate(lengths)}


def _revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _simulate(genome: dict, n: int, rng: random.Random, L=80, insert=250):
    """Exact substrings of the genome: ``n`` single-end reads, or mate pairs."""
    chroms = list(genome)
    se, m1, m2 = [], [], []
    for i in range(n):
        c = rng.choice(chroms)
        s = genome[c]
        p = rng.randrange(0, len(s) - insert)
        frag = s[p:p + insert]
        if rng.random() < 0.5:
            frag = _revcomp(frag)
        name = f"read{i}_{c}_{p}".encode()
        se.append((name, frag[:L].encode()))
        m1.append((name + b"/1", frag[:L].encode()))
        m2.append((name + b"/2", _revcomp(frag[-L:]).encode()))
    return se, m1, m2


@requires_pysam
@pytest.mark.skipif(not _TOOLS, reason="hisat2, hisat2-build and samtools needed")
def test_end_to_end_local_fastq_with_hisat2(tmp_path: Path):
    import subprocess
    import pysam
    rng = random.Random(12)
    genome = _random_genome(rng, [60_000, 40_000])
    gfa = tmp_path / "genome.fa"
    gfa.write_text("".join(f">{c}\n{s}\n" for c, s in genome.items()))
    idx_dir = tmp_path / "idx"
    idx_dir.mkdir()
    subprocess.run(["hisat2-build", "-q", str(gfa), str(idx_dir / "hisatidx")], check=True)

    se, m1, m2 = _simulate(genome, 2000, rng)
    _write_gz(tmp_path / "single.fq.gz", _fastq_bytes(se), members=3, rng=rng)
    _write_gz(tmp_path / "pair_1.fq.gz", _fastq_bytes(m1))
    _write_gz(tmp_path / "pair_2.fq.gz", _fastq_bytes(m2))

    out = tmp_path / "out"
    rc = cli.main(["run", "Test species", str(gfa), "--index", str(idx_dir / "hisatidx"),
                   "--outdir", str(out), "--threads", "2", "--keep-bam",
                   "--batch-size", "250", "--max-batches", "6", "--seed", "3",
                   "--tile-size", "1000",
                   "--fastq", str(tmp_path / "single.fq.gz"),
                   "--fastq", f"{tmp_path / 'pair_1.fq.gz'},{tmp_path / 'pair_2.fq.gz'}"])
    assert rc == 0
    bam = out / "VARUS.bam"
    assert bam.is_file()
    header, rows = read_manifest(out / MANIFEST_NAME)
    assert sum(r["n_batches"] for r in rows) == 6
    assert {r["accession"] for r in rows} <= {"single", "pair"}
    assert "pair" in {r["accession"] for r in rows}      # 2 mates > 1 read per spot
    assert header["local_fastq.pair"].endswith("pair_1.fq.gz,"
                                                + str((tmp_path / "pair_2.fq.gz").resolve()))
    # every read in the BAM comes from one of the sampled spot ranges
    with pysam.AlignmentFile(str(bam)) as af:
        n_aligned = 0
        for rd in af.fetch(until_eof=True):
            n_aligned += 1
            i = int(rd.query_name.split("_")[0][4:])
            acc = "pair" if rd.is_paired else "single"
            assert any(r["accession"] == acc and r["n"] <= i <= r["x"] for r in rows), rd.query_name
    assert n_aligned > 6 * 250 * 0.8
    assert (out / "RunStatistics.csv").is_file() and (out / "Coverage.csv").is_file()
    assert not (out / "batches").exists()               # batch FASTAs cleaned up

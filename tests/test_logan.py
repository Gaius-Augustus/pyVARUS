"""Tests for varus.logan (no network; HTTP is injected via ``opener``)."""

from __future__ import annotations

import io
import json
import math
import os
import random
import threading
import urllib.error
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pytest

from tests.conftest import requires_pysam, requires_zstandard
from varus import logan
from varus.logan import (
    LoganConfig,
    LoganContigs,
    LoganRunStats,
    apply_gate,
    check_availability,
    download_contigs,
    head_available,
    load_logan,
    load_run_stats,
    parse_logan_header,
    read_runlist,
    run_logan,
    sample_candidates,
    save_run_stats,
    scan_chunk_bam,
)
from varus.runlist import RunRecord, write_runlist


@pytest.fixture(autouse=True)
def _no_sleep(monkeypatch):
    """Backoff sleeps are pointless in tests."""
    monkeypatch.setattr(logan.time, "sleep", lambda s: None)


# ---------------------------------------------------------------------------
# parse_logan_header
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "line,name,ka",
    [
        (">SRR10620183_0 ka:f:431.586   L:-:10472:-  ", "SRR10620183_0", 431.586),
        (">SRR1_7 km:f:12.5", "SRR1_7", 12.5),
        (">SRR1_8 ka:f:3 L:+:1:- L:-:2:+\n", "SRR1_8", 3.0),
        ("SRR1_9 ka:f:2.0", "SRR1_9", 2.0),  # no leading '>'
        (">SRR1_3\n", "SRR1_3", math.nan),  # missing abundance
        (">SRR1_4 ka:f:abc", "SRR1_4", math.nan),  # unparseable
        (">SRR1_5 L:-:1:+", "SRR1_5", math.nan),  # only link fields
        (">   ", "", math.nan),
    ],
)
def test_parse_logan_header(line, name, ka):
    got_name, got_ka = parse_logan_header(line)
    assert got_name == name
    if math.isnan(ka):
        assert math.isnan(got_ka)
    else:
        assert got_ka == pytest.approx(ka)


# ---------------------------------------------------------------------------
# sample_candidates / read_runlist
# ---------------------------------------------------------------------------


def _rec(acc: str, bp: str = "") -> RunRecord:
    return RunRecord(accession=acc, total_spots=1000, total_bases=100_000,
                     avg_len=100.0, paired=False, colorspace=False,
                     platform="ILLUMINA", bioproject=bp)


def test_sample_candidates_deterministic_and_round_robin():
    recs = [_rec(f"SRR{i}", "PRJA") for i in range(5)] + \
           [_rec(f"ERR{i}", "PRJB") for i in range(3)] + \
           [_rec("DRR0", ""), _rec("DRR1", "")]
    out1 = sample_candidates(recs, 0)
    out2 = sample_candidates(list(reversed(recs)), 0)  # input order irrelevant
    assert out1 == out2
    assert len(out1) == len(recs)
    assert {r.accession for r in out1} == {r.accession for r in recs}
    # first round visits every bucket once: PRJA, PRJB, DRR0, DRR1 -> 4 buckets
    first_round = out1[:4]
    keys = [r.bioproject or r.accession for r in first_round]
    assert len(set(keys)) == 4
    # second round: only PRJA and PRJB have a second member
    second_round = out1[4:6]
    assert {r.bioproject for r in second_round} == {"PRJA", "PRJB"}
    # the remainder is all PRJA
    assert all(r.bioproject == "PRJA" for r in out1[8:])
    # max_candidates truncates the same ordering
    assert sample_candidates(recs, 3) == out1[:3]


def test_sample_candidates_orders_by_sha1():
    recs = [_rec("A", "P1"), _rec("B", "P1"), _rec("C", "P2")]
    out = sample_candidates(recs, 0)
    # within P1, order is by sha1(acc)
    p1 = [r.accession for r in out if r.bioproject == "P1"]
    assert p1 == sorted(p1, key=logan._sha1)


def test_read_runlist_matches_load_runs(tmp_path: Path):
    recs = [_rec("SRR1", "PRJA"), _rec("SRR2")]
    path = tmp_path / "Runlist.tsv"
    write_runlist(recs, path)
    with path.open("a") as f:
        f.write("SRR3\t10\t1000\t100.0\t0\t1\tABI_SOLID\tPRJC\n")  # colorspace
        f.write("bad\tline\n")
    got = read_runlist(path)
    assert [r.accession for r in got] == ["SRR1", "SRR2"]
    assert got[0].bioproject == "PRJA" and got[0].platform == "ILLUMINA"


# ---------------------------------------------------------------------------
# head_available / check_availability
# ---------------------------------------------------------------------------


class _Resp:
    def __init__(self, status: int, body: bytes = b""):
        self.status = status
        self._buf = io.BytesIO(body)
        self.closed = False

    def getcode(self):
        return self.status

    def read(self, n=-1):
        return self._buf.read(n)

    def close(self):
        self.closed = True


def _http_error(code: int):
    return urllib.error.HTTPError("http://x", code, "err", {}, None)


def _opener_from(script):
    """Return an opener that yields the scripted outcomes in order."""
    calls: List[str] = []
    it = iter(script)

    def opener(req, timeout):
        calls.append(req.full_url)
        assert req.get_method() in ("HEAD", "GET")
        item = next(it)
        if isinstance(item, Exception):
            raise item
        return item

    opener.calls = calls
    return opener


def test_head_available_200():
    op = _opener_from([_Resp(200)])
    assert head_available("SRR1", timeout=1, retries=3, opener=op) == ("available", 200)
    assert op.calls == [logan.LOGAN_CONTIGS_URL.format(acc="SRR1")]


def test_head_available_404_and_403_are_absent():
    op = _opener_from([_http_error(404)])
    assert head_available("SRR1", timeout=1, retries=3, opener=op) == ("absent", 404)
    op = _opener_from([_http_error(403)])
    assert head_available("SRR1", timeout=1, retries=3, opener=op) == ("absent", 403)
    assert len(op.calls) == 1  # no retry on absent


def test_head_available_500_retries_then_succeeds():
    op = _opener_from([_http_error(500), _http_error(503), _Resp(200)])
    assert head_available("SRR1", timeout=1, retries=3, opener=op) == ("available", 200)
    assert len(op.calls) == 3


def test_head_available_500_exhausts_retries():
    op = _opener_from([_http_error(500)] * 3)
    assert head_available("SRR1", timeout=1, retries=2, opener=op) == ("error", 500)
    assert len(op.calls) == 3  # retries + 1 attempts


def test_head_available_urlerror():
    op = _opener_from([urllib.error.URLError("dns"), urllib.error.URLError("dns")])
    assert head_available("SRR1", timeout=1, retries=1, opener=op) == ("error", 0)
    op = _opener_from([urllib.error.URLError("dns"), _Resp(200)])
    assert head_available("SRR1", timeout=1, retries=1, opener=op) == ("available", 200)


def test_head_available_custom_template():
    op = _opener_from([_Resp(200)])
    head_available("SRR1", timeout=1, retries=0, opener=op, url_template=logan.LOGAN_UNITIGS_URL)
    assert op.calls[0].endswith("/u/SRR1/SRR1.unitigs.fa.zst")


def test_check_availability_cache_resume(tmp_path: Path):
    cache = tmp_path / "availability.tsv"
    outcomes = {"SRR1": _Resp(200), "SRR2": _http_error(404), "SRR3": _http_error(500)}
    calls: List[str] = []

    def opener(req, timeout):
        acc = req.full_url.split("/")[-2]
        calls.append(acc)
        item = outcomes[acc]
        if isinstance(item, Exception):
            raise item
        return item

    res = check_availability(["SRR1", "SRR2", "SRR3"], cache, workers=2, timeout=1,
                             retries=0, opener=opener)
    assert res == {"SRR1": ("available", 200), "SRR2": ("absent", 404), "SRR3": ("error", 500)}
    assert cache.is_file()
    lines = cache.read_text().splitlines()
    assert lines[0] == "acc\tstatus\thttp\tchecked_at"
    assert len(lines) == 4
    # second call: only the "error" acc and the new acc are queried
    calls.clear()
    outcomes["SRR3"] = _Resp(200)
    outcomes["SRR4"] = _Resp(200)
    res = check_availability(["SRR1", "SRR2", "SRR3", "SRR4"], cache, workers=2,
                             timeout=1, retries=0, opener=opener)
    assert sorted(calls) == ["SRR3", "SRR4"]
    assert res["SRR1"] == ("available", 200)
    assert res["SRR3"] == ("available", 200)
    assert res["SRR4"] == ("available", 200)
    # third call: nothing to do
    calls.clear()
    check_availability(["SRR1", "SRR4"], cache, workers=2, timeout=1, retries=0, opener=opener)
    assert calls == []


# ---------------------------------------------------------------------------
# download_contigs
# ---------------------------------------------------------------------------


def _fasta_bytes(acc: str, n: int, ka_missing_at=None) -> bytes:
    out = []
    for i in range(n):
        if ka_missing_at is not None and i == ka_missing_at:
            out.append(f">{acc}_{i}\n")
        else:
            out.append(f">{acc}_{i} ka:f:{10.0 * (i + 1)}   L:-:1:-  \n")
        out.append("ACGT" * (i + 1) + "\n")
    return "".join(out).encode()


@requires_zstandard
def test_download_contigs_stream(tmp_path: Path, caplog):
    import zstandard

    raw = _fasta_bytes("SRR9", 3, ka_missing_at=1)
    comp = zstandard.ZstdCompressor().compress(raw)
    op = _opener_from([_Resp(200, comp)])
    with caplog.at_level("WARNING"):
        c = download_contigs("SRR9", tmp_path, timeout=1, retries=0, opener=op)
    assert isinstance(c, LoganContigs)
    assert c.acc == "SRR9"
    assert c.n_contigs == 3
    assert c.total_bp == 4 + 8 + 12
    assert c.bytes_downloaded == len(comp)
    assert c.ka.dtype == np.float32
    assert c.ka[0] == 10.0 and math.isnan(c.ka[1]) and c.ka[2] == 30.0
    assert sum("parseable abundance" in r.message for r in caplog.records) == 1
    text = c.fasta.read_text()
    assert text == ">SRR9_0\nACGT\n>SRR9_1\nACGTACGT\n>SRR9_2\nACGTACGTACGT\n"
    assert (tmp_path / "SRR9.ka.npy").is_file()
    assert not (tmp_path / "SRR9.contigs.fa.part").exists()

    # resume: no request is made
    op2 = _opener_from([])
    c2 = download_contigs("SRR9", tmp_path, timeout=1, retries=0, opener=op2)
    assert op2.calls == []
    assert c2.n_contigs == 3 and c2.total_bp == 24 and c2.bytes_downloaded == 0
    assert math.isnan(c2.ka[1])


@requires_zstandard
def test_download_contigs_retries_then_ok(tmp_path: Path):
    import zstandard

    comp = zstandard.ZstdCompressor().compress(_fasta_bytes("SRR8", 2))
    # first attempt: truncated stream -> decode error; second: 500; third ok
    op = _opener_from([_Resp(200, comp[:10]), _http_error(502), _Resp(200, comp)])
    c = download_contigs("SRR8", tmp_path, timeout=1, retries=2, opener=op)
    assert c.n_contigs == 2
    assert len(op.calls) == 3


def test_download_contigs_absent_raises(tmp_path: Path):
    op = _opener_from([_http_error(404)])
    with pytest.raises(logan.LoganAbsentError) as ei:
        download_contigs("SRR7", tmp_path, timeout=1, retries=3, opener=op)
    assert ei.value.http == 404
    assert len(op.calls) == 1


def test_download_contigs_gives_up(tmp_path: Path):
    op = _opener_from([_http_error(500)] * 2)
    with pytest.raises(logan.LoganDownloadError) as ei:
        download_contigs("SRR6", tmp_path, timeout=1, retries=1, opener=op)
    assert ei.value.http == 500
    assert not list(tmp_path.glob("SRR6*"))


# ---------------------------------------------------------------------------
# scan_chunk_bam
# ---------------------------------------------------------------------------


def _write_bam(path: Path, records, header=None):
    import pysam

    header = header or {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100_000},
                                                    {"SN": "chr2", "LN": 100_000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for qname, flag, ref_id, start, cigar in records:
            r = pysam.AlignedSegment(out.header)
            r.query_name = qname
            r.flag = flag
            qlen = sum(l for op, l in cigar if op in (0, 1, 4, 7, 8)) if cigar else 50
            r.query_sequence = "A" * qlen
            r.reference_id = ref_id
            r.reference_start = start
            r.mapping_quality = 60 if ref_id >= 0 else 0
            if cigar:
                r.cigartuples = cigar
            out.write(r)


@requires_pysam
@pytest.mark.parametrize("tile_weight", ["unit", "ka", "ka_len"])
def test_scan_chunk_bam(tmp_path: Path, tile_weight):
    bam = tmp_path / "chunk.bam"
    records = [
        # SRR1_0: spliced, tile 0 on chr1, aligned ref len 30+100+30=160
        ("SRR1_0", 0, 0, 0, [(0, 30), (3, 100), (0, 30)]),
        # SRR1_1: unspliced, tile 1 (POS 6001), 50M
        ("SRR1_1", 0, 0, 6_000, [(0, 50)]),
        # SRR1_1 secondary -> skipped
        ("SRR1_1", 256, 1, 100, [(0, 50)]),
        # SRR1_2 unmapped
        ("SRR1_2", 4, -1, -1, None),
        # SRR1_3 supplementary -> skipped
        ("SRR1_3", 2048, 1, 100, [(0, 50)]),
        # SRR2_0: spliced, chr2, tile 2 (POS 10001), aligned ref len 20+200+20=240 ; soft clip
        ("SRR2_0", 0, 1, 10_000, [(4, 5), (0, 20), (3, 200), (0, 20)]),
        # unknown run -> ignored with a warning
        ("XXX_0", 0, 0, 0, [(0, 50)]),
    ]
    _write_bam(bam, records)
    ka = {"SRR1": np.array([100.0, 5.0, math.nan, 1.0], dtype=np.float32),
          "SRR2": np.array([math.nan], dtype=np.float32),
          "SRR3": np.zeros(0, dtype=np.float32)}  # nothing aligned
    stats = scan_chunk_bam(bam, ka, {"SRR1": 4, "SRR2": 1, "SRR3": 0},
                           {"SRR1": 4000, "SRR2": 1000, "SRR3": 0},
                           tile_size=5000, ka_cap=50.0, tile_weight=tile_weight)
    assert set(stats) == {"SRR1", "SRR2", "SRR3"}
    s1, s2, s3 = stats["SRR1"], stats["SRR2"], stats["SRR3"]
    assert s1.n_contigs == 4 and s1.total_bp == 4000
    assert s1.n_aligned == 2 and s1.n_spliced == 1
    assert s1.aligned_bp == 60 + 50
    assert s1.mapped_pct == pytest.approx(50.0)
    assert s1.n_tiles == 2 and set(s1.tiles) == {("chr1", 0), ("chr1", 1)}
    assert s2.n_aligned == 1 and s2.n_spliced == 1 and s2.aligned_bp == 40
    assert set(s2.tiles) == {("chr2", 2)}
    assert s3.n_aligned == 0 and s3.n_tiles == 0 and s3.mapped_pct == 0.0

    w0, w1, w2 = s1.tiles[("chr1", 0)], s1.tiles[("chr1", 1)], s2.tiles[("chr2", 2)]
    if tile_weight == "unit":
        assert (w0, w1, w2) == (1.0, 1.0, 1.0)
    elif tile_weight == "ka":
        assert (w0, w1, w2) == (50.0, 5.0, 1.0)  # capped, plain, NaN->1
    else:
        assert w0 == pytest.approx(50.0 * 160 / 150)
        assert w1 == pytest.approx(5.0 * 1.0)  # 50 < 150 -> factor 1
        assert w2 == pytest.approx(1.0 * 240 / 150)
    assert s1.tile_mass == pytest.approx(w0 + w1)
    # introns: 1-based inclusive, weight = min(ka, cap), NaN -> 1
    assert s1.introns == {("chr1", 31, 130, "."): 50.0}
    assert s2.introns == {("chr2", 10_021, 10_220, "."): 1.0}


@requires_pysam
def test_scan_chunk_bam_tile_ka_cap(tmp_path: Path):
    """tile_ka_cap changes the tile weights only; introns keep ka_cap."""
    bam = tmp_path / "chunk.bam"
    _write_bam(bam, [("SRR1_0", 0, 0, 0, [(0, 30), (3, 100), (0, 30)])])
    ka = {"SRR1": np.array([400.0], dtype=np.float32)}
    args = (bam, ka, {"SRR1": 1}, {"SRR1": 160})
    for cap, want in ((None, 50.0), (0.0, 400.0), (100.0, 100.0)):
        st = scan_chunk_bam(*args, tile_size=5000, ka_cap=50.0, tile_weight="ka",
                            tile_ka_cap=cap)["SRR1"]
        assert st.tiles[("chr1", 0)] == pytest.approx(want)
        assert st.introns == {("chr1", 31, 130, "."): 50.0}


@requires_pysam
def test_scan_chunk_bam_divergence(tmp_path: Path):
    import pysam

    bam = tmp_path / "chunk.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100_000}]}
    # (qname, aligned length, de) ; SRR2_0 has only NM (fallback NM / len)
    recs = [("SRR1_0", 100, 0.20), ("SRR1_1", 300, 0.01), ("SRR1_2", 50, 0.30),
            ("SRR2_0", 200, None)]
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for i, (q, n, de) in enumerate(recs):
            r = pysam.AlignedSegment(out.header)
            r.query_name, r.flag, r.reference_id = q, 0, 0
            r.reference_start, r.mapping_quality = 1000 * i, 60
            r.query_sequence = "A" * n
            r.cigartuples = [(0, n)]
            if de is None:
                r.set_tag("NM", 10)
            else:
                r.set_tag("de", de, value_type="f")
            out.write(r)
    ka = {a: np.ones(3, dtype=np.float32) for a in ("SRR1", "SRR2", "SRR3")}
    stats = scan_chunk_bam(bam, ka, {}, {}, tile_size=5000, ka_cap=50.0, tile_weight="unit")
    # weighted median: 300 of 450 bp at 0.01
    assert stats["SRR1"].divergence == pytest.approx(0.01)
    assert stats["SRR2"].divergence == pytest.approx(10 / 200)
    assert math.isnan(stats["SRR3"].divergence)


@requires_pysam
def test_scan_chunk_bam_rejects_bad_args(tmp_path: Path):
    with pytest.raises(ValueError):
        scan_chunk_bam(tmp_path / "x.bam", {}, {}, {}, tile_size=5000, ka_cap=50, tile_weight="foo")
    with pytest.raises(ValueError):
        scan_chunk_bam(tmp_path / "x.bam", {}, {}, {}, tile_size=0, ka_cap=50, tile_weight="unit")


# ---------------------------------------------------------------------------
# apply_gate
# ---------------------------------------------------------------------------


def _st(acc, n_contigs, n_tiles, status="pending"):
    return LoganRunStats(acc=acc, n_contigs=n_contigs, n_tiles=n_tiles,
                         tiles={("chr1", i): 1.0 for i in range(n_tiles)}, status=status)


def test_apply_gate():
    stats = {
        "GOOD": _st("GOOD", 10_000, 2400),
        "OK": _st("OK", 5_000, 300),
        "FOREIGN": _st("FOREIGN", 8_000, 150),
        "FEW": _st("FEW", 50, 5000),  # too few contigs, must not set the max
        "ZERO": _st("ZERO", 500, 0),
        "ABSENT": LoganRunStats(acc="ABSENT", status="absent"),
        "ERR": LoganRunStats(acc="ERR", status="error"),
    }
    apply_gate(stats, min_contigs=100, min_tiles_frac=0.10)
    assert stats["GOOD"].status == "accepted"
    assert stats["OK"].status == "accepted"     # 300 >= 240
    assert stats["FOREIGN"].status == "rejected"  # 150 < 240
    assert stats["FEW"].status == "too_few_contigs"
    assert stats["ZERO"].status == "rejected"
    assert stats["ABSENT"].status == "absent"
    assert stats["ERR"].status == "error"


def test_apply_gate_divergence():
    # another species: broad tile coverage, but 15 % divergent contigs
    other = _st("OTHER", 30_000, 5000)
    other.divergence = 0.15
    good = _st("GOOD", 10_000, 2400)
    good.divergence = 0.0
    unknown = _st("UNKNOWN", 10_000, 2000)  # NaN (old cache) passes
    stats = {"OTHER": other, "GOOD": good, "UNKNOWN": unknown}
    apply_gate(stats, min_contigs=100, min_tiles_frac=0.9)
    assert stats["OTHER"].status == "rejected"
    # OTHER must not set max_tiles: 0.9 * 2400 = 2160 -> GOOD in, UNKNOWN out
    assert stats["GOOD"].status == "accepted"
    assert stats["UNKNOWN"].status == "rejected"
    # 0 switches the divergence check off
    apply_gate(stats, min_contigs=100, min_tiles_frac=0.1, max_divergence=0)
    assert stats["OTHER"].status == "accepted"


def test_apply_gate_no_eligible_runs():
    stats = {"FEW": _st("FEW", 10, 100)}
    apply_gate(stats, min_contigs=100, min_tiles_frac=0.1)
    assert stats["FEW"].status == "too_few_contigs"


# ---------------------------------------------------------------------------
# save / load run stats
# ---------------------------------------------------------------------------


def test_save_load_run_stats_roundtrip(tmp_path: Path):
    st = LoganRunStats(
        acc="SRR1", n_contigs=10, total_bp=12345, n_aligned=7, aligned_bp=999,
        mapped_pct=70.0, n_spliced=3, n_tiles=2, tile_mass=3.5,
        tiles={("chr1", 0): 1.5, ("chr_2", 12): 2.0},
        introns={("chr1", 31, 130, "."): 50.0, ("chr_2", 5, 10, "."): 1.0},
        status="accepted", rank=2, gain=12.5, http=200, bioproject="PRJX",
    )
    jpath, npath = save_run_stats(st, tmp_path / "runs")
    assert jpath.is_file() and npath.is_file()
    back = load_run_stats("SRR1", tmp_path / "runs")
    assert back == st
    assert load_run_stats("NOPE", tmp_path / "runs") is None
    # empty tiles/introns also round-trip
    empty = LoganRunStats(acc="SRR2", status="rejected")
    save_run_stats(empty, tmp_path / "runs")
    assert load_run_stats("SRR2", tmp_path / "runs") == empty
    st.divergence = 0.0123
    save_run_stats(st, tmp_path / "runs")
    assert load_run_stats("SRR1", tmp_path / "runs").divergence == pytest.approx(0.0123)


# ---------------------------------------------------------------------------
# end-to-end run_logan with the network and aligner stubbed out
# ---------------------------------------------------------------------------

TILE = 5000
N_TILES = 60
GOOD = {"SRR1": 400, "SRR2": 400, "ERR3": 300}  # acc -> n_contigs (same bioproject SRR1/SRR2)
FOREIGN = {"DRR4": 500}
TOO_FEW = {"SRR5": 20}
EMPTY = {"SRR7": 200}  # contigs but nothing aligns
ABSENT = ["SRR6"]
# Tiles covered by the GOOD runs: SRR2 all 60, ERR3 50, SRR1 35.
COVER = {"SRR1": 35, "SRR2": N_TILES, "ERR3": 50}


def _make_genome(path: Path) -> None:
    rng = random.Random(3)
    seq = bytearray(rng.choice(b"ACGT") for _ in range(N_TILES * TILE))
    for t in range(N_TILES):
        s = t * TILE + 10
        seq[s + 30:s + 32] = b"GT"
        seq[s + 228:s + 230] = b"AG"
    with path.open("w") as f:
        f.write(">chr1\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i + 80].decode() + "\n")


def _fake_download_factory(counts: Dict[str, int], calls: List[str]):
    def fake_download(acc, dest_dir, *, timeout, retries, url_template=None, opener=None):
        calls.append(acc)
        if acc in ABSENT:
            raise logan.LoganAbsentError("gone", http=404)
        n = counts[acc]
        dest_dir = Path(dest_dir)
        dest_dir.mkdir(parents=True, exist_ok=True)
        fasta = dest_dir / f"{acc}.contigs.fa"
        with fasta.open("w") as f:
            for i in range(n):
                f.write(f">{acc}_{i}\n" + "ACGT" * 50 + "\n")
        ka = np.full(n, 20.0, dtype=np.float32)
        np.save(dest_dir / f"{acc}.ka.npy", ka)
        return LoganContigs(acc=acc, fasta=fasta, ka=ka, n_contigs=n, total_bp=200 * n,
                            bytes_downloaded=1_000_000, seconds=0.5)
    return fake_download


def _fake_align(queries, *, index, out_bam, threads=4, max_intron=20000, **kw):
    """Emit a BAM whose alignments depend on the run's role."""
    records = []
    for q in queries:
        acc = Path(q).name.split(".")[0]
        names = [l[1:].strip() for l in Path(q).read_text().splitlines() if l.startswith(">")]
        for name in names:
            i = int(name.rsplit("_", 1)[1])
            if acc in GOOD:
                t = i % N_TILES
                if t >= COVER[acc]:
                    continue
                s = t * TILE + 10
                if (i // N_TILES) % 2 == 0:  # every tile gets spliced contigs
                    records.append((name, 0, 0, s, [(0, 30), (3, 200), (0, 30)]))
                else:
                    records.append((name, 0, 0, s, [(0, 60)]))
            elif acc in FOREIGN:
                if i < 3:
                    records.append((name, 0, 0, (i % 2) * TILE + 10, [(0, 60)]))
                else:
                    records.append((name, 4, -1, -1, None))
            elif acc in TOO_FEW:
                records.append((name, 0, 0, 10, [(0, 60)]))
    Path(out_bam).parent.mkdir(parents=True, exist_ok=True)
    _write_bam(Path(out_bam), records,
               header={"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": N_TILES * TILE}]})
    return Path(out_bam)


class _FakeAligner:
    """Stand-in for the minimap2 Popen: writes a prepared BAM into the pipe."""

    def __init__(self, fd: int, data: bytes) -> None:
        self._t = threading.Thread(target=self._run, args=(os.dup(fd), data))
        self._t.start()

    @staticmethod
    def _run(fd: int, data: bytes) -> None:
        with os.fdopen(fd, "wb") as f:
            f.write(data)

    def wait(self) -> int:
        self._t.join()
        return 0

    def kill(self) -> None:
        pass


def _fake_start_alignment(queries, *, index, stdout, threads=4, max_intron=20000, log_path, **kw):
    """Streamed variant of :func:`_fake_align`: same BAM, sent through the pipe."""
    tmp = Path(log_path).with_suffix(".fake.bam")
    _fake_align(queries, index=index, out_bam=tmp, threads=threads, max_intron=max_intron)
    data = tmp.read_bytes()
    tmp.unlink()
    Path(log_path).write_text("[M::main] Real time: 0.010 sec; CPU: 0.020 sec; Peak RSS: 0.001 GB\n")
    return _FakeAligner(stdout, data)


def _e2e_setup(tmp_path: Path, monkeypatch, head_status="available"):
    genome = tmp_path / "genome.fa"
    _make_genome(genome)
    recs = [_rec("SRR1", "PRJA"), _rec("SRR2", "PRJA"), _rec("ERR3", "PRJB"),
            _rec("DRR4", "PRJC"), _rec("SRR5", "PRJD"), _rec("SRR6", "PRJE"),
            _rec("SRR7", "PRJF")]
    runlist = tmp_path / "Runlist.tsv"
    write_runlist(recs, runlist)
    counts = {**GOOD, **FOREIGN, **TOO_FEW, **EMPTY}
    dl_calls: List[str] = []
    head_calls: List[str] = []

    def fake_head(acc, *, timeout, retries, url_template=None, opener=None):
        head_calls.append(acc)
        if head_status != "available":
            return head_status, 500 if head_status == "error" else 404
        return ("absent", 404) if acc in ABSENT else ("available", 200)

    monkeypatch.setattr(logan, "head_available", fake_head)
    monkeypatch.setattr(logan, "download_contigs", _fake_download_factory(counts, dl_calls))
    monkeypatch.setattr(logan, "align_contigs_minimap2", _fake_align)
    monkeypatch.setattr(logan, "start_contig_alignment", _fake_start_alignment)
    monkeypatch.setattr(logan, "build_minimap2_index",
                        lambda g, o, t: (Path(o).mkdir(parents=True, exist_ok=True),
                                         Path(o) / "mm2idx.mmi")[1])
    cfg = LoganConfig(genome=genome, runlist=runlist, outdir=tmp_path / "out",
                      threads=1, download_workers=2, max_candidates=6, chunk_runs=2,
                      min_contigs=100, min_tiles_frac=0.10, select_top=2,
                      batch_size=1000, prior_batches=1.0, tile_size=TILE)
    return cfg, recs, dl_calls, head_calls


@requires_pysam
def test_run_logan_end_to_end(tmp_path: Path, monkeypatch):
    cfg, recs, dl_calls, head_calls = _e2e_setup(tmp_path, monkeypatch)
    rc = run_logan(cfg)
    assert rc == 0
    ldir = cfg.outdir / "logan"

    # availability: SRR7 is beyond max_candidates=6 unless topped up; SRR6 is
    # absent so the top-up pulls SRR7 in as the 6th available candidate.
    assert (ldir / "availability.tsv").is_file()
    assert set(head_calls) == {r.accession for r in recs}
    assert sorted(dl_calls) == sorted(["SRR1", "SRR2", "ERR3", "DRR4", "SRR5", "SRR7"])

    # ranking table
    ranking = ldir / "LoganRanking.tsv"
    rows = [l.split("\t") for l in ranking.read_text().splitlines()]
    header = rows[0]
    assert header == ["acc", "bioproject", "status", "http", "n_contigs", "contig_mb",
                      "mapped_pct", "yield_pct", "divergence", "n_tiles", "tile_mass", "n_spliced",
                      "n_introns", "selected_rank", "gain", "cumulative_S"]
    by_acc = {r[0]: dict(zip(header, r)) for r in rows[1:]}
    assert set(by_acc) == {r.accession for r in recs}
    assert by_acc["SRR1"]["status"] == "accepted"
    assert by_acc["SRR2"]["status"] == "accepted"
    assert by_acc["ERR3"]["status"] == "accepted"
    assert by_acc["DRR4"]["status"] == "rejected"
    assert by_acc["SRR5"]["status"] == "too_few_contigs"
    assert by_acc["SRR6"]["status"] == "absent" and by_acc["SRR6"]["http"] == "404"
    assert by_acc["SRR6"]["n_contigs"] == ""
    assert int(by_acc["SRR2"]["n_tiles"]) == 60 and int(by_acc["SRR1"]["n_tiles"]) == 35
    assert int(by_acc["ERR3"]["n_tiles"]) == 50
    assert int(by_acc["DRR4"]["n_tiles"]) == 2
    assert float(by_acc["DRR4"]["mapped_pct"]) == pytest.approx(100 * 3 / 500)
    # SRR7 has contigs but none aligns -> 0 tiles -> rejected
    assert by_acc["SRR7"]["status"] == "rejected"
    assert int(by_acc["SRR7"]["n_tiles"]) == 0
    # selection: SRR2 covers all 60 tiles -> rank 1; ERR3 (50 tiles) beats SRR1 (35)
    assert by_acc["SRR2"]["selected_rank"] == "1"
    assert by_acc["ERR3"]["selected_rank"] == "2"
    assert by_acc["SRR1"]["selected_rank"] == ""
    assert float(by_acc["ERR3"]["cumulative_S"]) > float(by_acc["SRR2"]["cumulative_S"])

    # Runlist.logan.tsv: ranked, then accepted by n_tiles desc, then unprocessed
    out_rl = cfg.outdir / "Runlist.logan.tsv"
    from varus.controller import load_runs
    states = load_runs(out_rl, 1000, random.Random(0))
    accs = [s.record.accession for s in states]
    assert accs == ["SRR2", "ERR3", "SRR1", "SRR5", "SRR6"]
    assert not out_rl.read_text().lstrip().startswith("#")

    # introns / splice sites / junc bed
    gff = (ldir / "logan_introns.gff").read_text().splitlines()
    assert gff and all("\tintron\t" in l for l in gff)
    assert all(l.split("\t")[6] == "+" for l in gff)
    starts = {int(l.split("\t")[3]) for l in gff}
    assert starts == {t * TILE + 10 + 31 for t in range(N_TILES)}
    mults = [int(l.split("\t")[5]) for l in gff]
    assert min(mults) >= 20  # ka weight 20 per contig, summed over runs
    ss = (ldir / "logan.splice_sites").read_text().splitlines()
    assert len(ss) == len(gff)
    bed = (ldir / "logan.junc.bed").read_text().splitlines()
    assert len(bed) == len(gff) and len(bed[0].split("\t")) == 12

    # tiles table
    import gzip
    with gzip.open(ldir / "logan_tiles.tsv.gz", "rt") as f:
        tl = f.read().splitlines()
    assert tl[0] == "acc\tchrom\ttile_idx\tweight"
    tile_accs = {l.split("\t")[0] for l in tl[1:]}
    assert tile_accs == {"SRR1", "SRR2", "ERR3"}
    assert sum(1 for l in tl[1:] if l.startswith("SRR2\t")) == 60

    # summary
    summ = json.loads((ldir / "logan_summary.json").read_text())
    assert summ["params"]["tile_size"] == TILE
    assert summ["params"]["genome"] == str(cfg.genome)
    assert summ["counts"] == {"accepted": 3, "rejected": 2, "too_few_contigs": 1, "absent": 1}
    assert [c["acc"] for c in summ["coverage_curve"]] == ["SRR2", "ERR3"]
    assert summ["coverage_curve"][0]["cumulative_tiles"] == 60
    assert summ["download_mb"] == pytest.approx(6.0)
    for k in ("head", "align", "scan", "select", "total"):
        assert k in summ["timings_s"]

    # contigs deleted, no BAM kept
    assert list((ldir / "contigs").glob("*.fa")) == []
    assert not (ldir / "bams").exists()
    assert not (ldir / "LOGAN.bam").exists()
    # per-run stats present for every processed run
    assert sorted(p.stem for p in (ldir / "runs").glob("*.json")) == \
        sorted(["SRR1", "SRR2", "ERR3", "DRR4", "SRR5", "SRR7"])

    # load_logan round trip
    prior = load_logan(ldir)
    assert prior.status["SRR2"] == "accepted" and prior.status["SRR6"] == "absent"
    assert prior.rank == {"SRR1": 0, "SRR2": 1, "ERR3": 2, "DRR4": 0, "SRR5": 0, "SRR6": 0, "SRR7": 0}
    assert set(prior.tiles) == {"SRR1", "SRR2", "ERR3"}
    assert len(prior.tiles["SRR2"]) == 60
    assert prior.introns == ldir / "logan_introns.gff"
    assert prior.splice_sites == ldir / "logan.splice_sites"
    assert prior.junc_bed == ldir / "logan.junc.bed"
    assert prior.params["select_top"] == 2

    # resume: a second run does not download again but produces the same ranking
    dl_calls.clear()
    rc = run_logan(cfg)
    assert rc == 0
    assert dl_calls == []
    assert ranking.read_text() == "\n".join("\t".join(r) for r in rows) + "\n"


@requires_pysam
def test_run_logan_keep_contigs_and_bam(tmp_path: Path, monkeypatch):
    cfg, *_ = _e2e_setup(tmp_path, monkeypatch)
    cfg.keep_contigs = True
    cfg.write_bam = True
    merged = {}
    monkeypatch.setattr(logan, "_samtools_merge",
                        lambda bams, out, threads, samtools="samtools": merged.setdefault("bams", list(bams)))
    assert run_logan(cfg) == 0
    ldir = cfg.outdir / "logan"
    assert len(list((ldir / "contigs").glob("*.contigs.fa"))) == 6
    chunks = sorted((ldir / "bams").glob("chunk_*.bam"))
    assert len(chunks) == 3  # 6 runs / chunk_runs=2
    assert merged["bams"] == chunks


def test_run_logan_all_absent_returns_3(tmp_path: Path, monkeypatch):
    cfg, *_ = _e2e_setup(tmp_path, monkeypatch, head_status="absent")
    assert run_logan(cfg) == 3
    assert (cfg.outdir / "Runlist.logan.tsv").is_file()


def test_run_logan_unreachable_returns_4(tmp_path: Path, monkeypatch):
    cfg, *_ = _e2e_setup(tmp_path, monkeypatch, head_status="error")
    assert run_logan(cfg) == 4


@requires_pysam
def test_run_logan_parallel_scan_matches_serial(tmp_path: Path, monkeypatch):
    """Scanning chunk BAMs in worker processes must not change any output."""
    outputs = {}
    for workers in (0, 2):
        sub = tmp_path / f"w{workers}"
        sub.mkdir()
        cfg, *_ = _e2e_setup(sub, monkeypatch)
        cfg.scan_workers = workers
        assert run_logan(cfg) == 0
        ldir = cfg.outdir / "logan"
        outputs[workers] = (
            (ldir / "LoganRanking.tsv").read_text(),
            (ldir / "logan_introns.gff").read_text(),
            (cfg.outdir / "Runlist.logan.tsv").read_text(),
        )
        assert not (ldir / "tmp_bams").exists()  # chunk BAMs reaped after the scan
    assert outputs[0] == outputs[2]


def test_load_logan_tolerates_missing_files(tmp_path: Path):
    prior = load_logan(tmp_path / "nonexistent")
    assert prior.status == {} and prior.rank == {} and prior.tiles == {}
    assert prior.introns is None and prior.params == {}


@requires_pysam
def test_run_logan_minimap2_threads_leave_room_for_scanners(tmp_path: Path, monkeypatch):
    cfg, *_ = _e2e_setup(tmp_path, monkeypatch)
    cfg.threads = 16                       # scan_workers 2, align_groups 3 -> 3 scanners
    cfg.align_groups = 3                   # auto would give 1 at 16 threads
    cfg.chunk_runs = 3
    seen: List[Tuple[int, int]] = []

    def rec_start(queries, **kw):
        seen.append((len(queries), kw["threads"]))
        return _fake_start_alignment(queries, **kw)

    monkeypatch.setattr(logan, "start_contig_alignment", rec_start)
    assert run_logan(cfg) == 0
    # 16 threads - 3 scanners - 1 (main + downloads), capped at a quarter -> 12,
    # split over the minimap2 groups of each chunk (one run per group here)
    assert seen and all(q == 1 for q, _ in seen)
    assert {t for _, t in seen} <= {12 // 3, 12 // 2, 12}
    assert 12 // 3 in {t for _, t in seen}


@requires_pysam
def test_run_logan_align_groups_give_identical_results(tmp_path: Path, monkeypatch):
    out = {}
    for g in (1, 3):
        (tmp_path / f"g{g}").mkdir()
        cfg, *_ = _e2e_setup(tmp_path / f"g{g}", monkeypatch)
        cfg.align_groups = g
        cfg.chunk_runs = 6
        assert run_logan(cfg) == 0
        out[g] = (cfg.outdir / "logan" / "LoganRanking.tsv").read_text()
    assert out[1] == out[3]


def test_auto_align_groups_scale_with_threads():
    from varus.logan import auto_align_groups
    assert [auto_align_groups(t) for t in (1, 4, 8, 16, 26, 27, 32, 48, 64, 256)] == \
        [1, 1, 1, 1, 1, 2, 2, 3, 4, 4]


def test_max_groups_for_memory(tmp_path: Path):
    from varus.logan import available_memory_bytes, max_groups_for_memory
    mmi = tmp_path / "g.mmi"
    with open(mmi, "wb") as fh:
        fh.truncate(8 * 2**30)                       # sparse 8 GiB "index"
    gib = 2**30
    assert max_groups_for_memory(mmi, avail=100 * gib) == 6   # 60 GiB / 10 GiB
    assert max_groups_for_memory(mmi, avail=30 * gib) == 1
    assert max_groups_for_memory(mmi, avail=1 * gib) == 1     # never 0
    avail = available_memory_bytes()
    assert avail is None or avail > 0

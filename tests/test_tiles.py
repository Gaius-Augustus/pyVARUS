"""Tests for varus.tiles. BAM-touching tests require pysam."""

from __future__ import annotations

from pathlib import Path

import pytest

from varus import tiles
from tests.conftest import requires_pysam


def test_count_umrs_rejects_zero_tile_size(tmp_path: Path):
    with pytest.raises(ValueError, match="tile_size"):
        tiles.count_umrs_per_tile(tmp_path / "x.bam", tile_size=0)


def test_count_bam_stats_rejects_zero_tile_size(tmp_path: Path):
    with pytest.raises(ValueError, match="tile_size"):
        tiles.count_bam_stats(tmp_path / "x.bam", tile_size=0)


@requires_pysam
def test_count_umrs_basic(tmp_path: Path):
    """One single-mapper, one multi-mapper, one paired-end pair on same tile."""
    import pysam

    bam_path = tmp_path / "t.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 100_000}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        # r1: single alignment in tile 0 (POS=1) -> UMR, +1 to ('chr1', 0)
        r1 = pysam.AlignedSegment(out.header)
        r1.query_name = "r1"
        r1.query_sequence = "A" * 50
        r1.flag = 0
        r1.reference_id = 0
        r1.reference_start = 0  # SAM POS=1, tile = 1//5000 = 0
        r1.mapping_quality = 60
        r1.cigartuples = [(0, 50)]
        r1.query_qualities = pysam.qualitystring_to_array("I" * 50)
        out.write(r1)

        # r2: two alignments to different tiles -> not UMR
        for ref_start in (10, 6_000):  # tile 0 and tile 1
            r = pysam.AlignedSegment(out.header)
            r.query_name = "r2"
            r.query_sequence = "C" * 50
            r.flag = 0 if ref_start == 10 else 256  # secondary
            r.reference_id = 0
            r.reference_start = ref_start
            r.mapping_quality = 60
            r.cigartuples = [(0, 50)]
            r.query_qualities = pysam.qualitystring_to_array("I" * 50)
            out.write(r)

        # r3: paired-end pair both on tile 2 (POS in [10001, 15000]) -> UMR, +1
        for flag, ref_start in [(99, 10_010), (147, 10_200)]:
            r = pysam.AlignedSegment(out.header)
            r.query_name = "r3"
            r.query_sequence = "G" * 50
            r.flag = flag
            r.reference_id = 0
            r.reference_start = ref_start
            r.mapping_quality = 60
            r.cigartuples = [(0, 50)]
            r.query_qualities = pysam.qualitystring_to_array("I" * 50)
            out.write(r)

    counts = tiles.count_umrs_per_tile(bam_path, tile_size=5000)
    assert counts.get(("chr1", 0)) == 1   # r1
    assert counts.get(("chr1", 2)) == 1   # r3 pair
    # r2 split across tiles 0 and 1 -> not counted
    assert ("chr1", 1) not in counts
    # 3 distinct reads, but only 2 are UMR
    assert sum(counts.values()) == 2


@requires_pysam
def test_count_bam_stats_spliced_count(tmp_path: Path):
    """count_bam_stats tracks spliced reads (N-op in CIGAR)."""
    import pysam

    bam_path = tmp_path / "s.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100_000}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        # unspliced read
        r1 = pysam.AlignedSegment(out.header)
        r1.query_name = "r1"
        r1.query_sequence = "A" * 50
        r1.flag = 0
        r1.reference_id = 0
        r1.reference_start = 0
        r1.mapping_quality = 60
        r1.cigartuples = [(0, 50)]   # 50M – no intron
        r1.query_qualities = pysam.qualitystring_to_array("I" * 50)
        out.write(r1)

        # spliced read with an N-op
        r2 = pysam.AlignedSegment(out.header)
        r2.query_name = "r2"
        r2.query_sequence = "G" * 60
        r2.flag = 0
        r2.reference_id = 0
        r2.reference_start = 0
        r2.mapping_quality = 60
        r2.cigartuples = [(0, 30), (3, 100), (0, 30)]  # 30M 100N 30M
        r2.query_qualities = pysam.qualitystring_to_array("I" * 60)
        out.write(r2)

    stats = tiles.count_bam_stats(bam_path, tile_size=5000)
    assert stats.n_reads == 2
    assert stats.n_spliced == 1   # only r2 has an N-op


@requires_pysam
def test_count_umrs_three_hits_same_tile_excluded(tmp_path: Path):
    """RNAread::UMR rejects reads with > 2 hits to a single tile."""
    import pysam

    bam_path = tmp_path / "t.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100_000}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for i, ref_start in enumerate([0, 100, 200]):
            r = pysam.AlignedSegment(out.header)
            r.query_name = "rN"
            r.query_sequence = "A" * 50
            r.flag = 0 if i == 0 else 256
            r.reference_id = 0
            r.reference_start = ref_start
            r.mapping_quality = 60
            r.cigartuples = [(0, 50)]
            r.query_qualities = pysam.qualitystring_to_array("I" * 50)
            out.write(r)

    counts = tiles.count_umrs_per_tile(bam_path, tile_size=5000)
    assert counts == {}


@requires_pysam
def test_scan_batch_bam_matches_two_passes(tmp_path: Path):
    """scan_batch_bam == count_bam_stats + extract_introns_from_bam."""
    import pysam

    from varus.introns import extract_introns_from_bam

    bam_path = tmp_path / "m.bam"
    header = {"HD": {"VN": "1.6"},
              "SQ": [{"SN": "chr1", "LN": 100_000}, {"SN": "chr2", "LN": 100_000}]}
    recs = [
        # (name, flag, ref, start, cigar)
        ("a", 99, 0, 100, [(0, 30), (3, 200), (0, 30)]),     # pair, spliced
        ("a", 147, 0, 400, [(0, 60)]),
        ("b", 0, 0, 9_000, [(0, 20), (3, 50), (0, 20), (3, 70), (0, 20)]),
        ("b", 256, 1, 500, [(0, 20), (3, 50), (0, 40)]),     # secondary, other tile
        ("c", 0, 0, 20_000, [(0, 30), (3, 200), (0, 30)]),
        ("d", 0, 1, 100, [(0, 30), (3, 200), (0, 30)]),
        ("d", 2048, 1, 100, [(0, 30), (3, 200), (0, 30)]),   # supplementary, same tile
        ("e", 4, -1, -1, None),                              # unmapped
        ("f", 0, 0, 100, [(0, 30), (3, 200), (0, 30)]),      # repeats a's intron
    ]
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for name, flag, ref, start, cigar in recs:
            r = pysam.AlignedSegment(out.header)
            r.query_name = name
            r.flag = flag
            r.reference_id = ref
            r.reference_start = start
            r.mapping_quality = 60 if ref >= 0 else 0
            qlen = sum(l for op, l in cigar if op in (0, 1, 4)) if cigar else 50
            r.query_sequence = "A" * qlen
            if cigar:
                r.cigartuples = cigar
            out.write(r)

    stats, introns = tiles.scan_batch_bam(bam_path, tile_size=5000)
    assert introns.counts == extract_introns_from_bam(bam_path).counts
    assert introns.counts[("chr1", 131, 330, ".")] == 2       # a + f
    assert introns.counts[("chr2", 131, 330, ".")] == 2       # d primary + supplementary
    assert introns.counts[("chr2", 521, 570, ".")] == 1       # b's secondary
    assert stats.n_reads == 5                                 # e is unmapped
    assert stats.n_spliced == 5                               # a, b, c, d, f
    assert stats.umr_counts == {("chr1", 0): 2, ("chr1", 4): 1, ("chr2", 0): 1}


def _random_hisat_bam(path: Path, seed: int = 0, n_reads: int = 3000) -> None:
    """Coordinate-sorted BAM with pairs, split pairs, multimappers, NH tags."""
    import random

    import pysam

    rng = random.Random(seed)
    refs = [("chr1", 230_000), ("chr2", 120_000), ("chrM", 9_000)]
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": n, "LN": l} for n, l in refs]}
    recs = []

    def seg(name, flag, ref, pos, nh, spliced, mate=None):
        recs.append((name, flag, ref, pos, nh, spliced, mate))

    for i in range(n_reads):
        name = f"r{i}"
        kind = rng.random()
        ref = rng.randrange(len(refs))
        pos = rng.randrange(0, refs[ref][1] - 500)
        if rng.random() < 0.05:   # start exactly on a 5 kb boundary
            pos = (pos // 5000) * 5000
        spl = rng.random() < 0.3
        if kind < 0.35:                                   # single end, unique
            seg(name, 0, ref, pos, 1, spl)
        elif kind < 0.7:                                  # proper pair, nearby
            p2 = min(refs[ref][1] - 200, pos + rng.randrange(50, 12_000))
            seg(name, 99, ref, pos, 1, spl, (ref, p2))
            seg(name, 147, ref, p2, 1, False, (ref, pos))
        elif kind < 0.8:                                  # pair across chromosomes
            ref2 = (ref + 1) % len(refs)
            p2 = rng.randrange(0, refs[ref2][1] - 500)
            seg(name, 97, ref, pos, 1, spl, (ref2, p2))
            seg(name, 145, ref2, p2, 1, False, (ref, pos))
        elif kind < 0.85:                                 # mate unmapped
            seg(name, 73, ref, pos, 1, spl)
        elif kind < 0.93:                                 # multimapper, 2 loci
            p2 = rng.choice([pos + rng.randrange(0, 300), rng.randrange(0, refs[ref][1] - 500)])
            seg(name, 0, ref, pos, 2, spl)
            seg(name, 256, ref, min(p2, refs[ref][1] - 500), 2, spl)
        else:                                             # no NH tag
            seg(name, 0, ref, pos, None, spl)

    unsorted = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for name, flag, ref, pos, nh, spliced, mate in recs:
            r = pysam.AlignedSegment(out.header)
            r.query_name, r.flag, r.reference_id, r.reference_start = name, flag, ref, pos
            r.mapping_quality = 60
            cig = [(0, 30), (3, 150), (0, 30)] if spliced else [(0, 60)]
            r.cigartuples = cig
            r.query_sequence = "A" * sum(l for op, l in cig if op == 0)
            if mate is not None:
                r.next_reference_id, r.next_reference_start = mate
            if nh is not None:
                r.set_tag("NH", nh)
            out.write(r)
    pysam.sort("-o", str(path), str(unsorted))
    unsorted.unlink()


@requires_pysam
@pytest.mark.parametrize("n_parts", [1, 2, 3, 7, 16, 64])
def test_scan_batch_bam_parallel_equals_one_pass(tmp_path: Path, n_parts):
    from concurrent.futures import ThreadPoolExecutor

    bam = tmp_path / "b.bam"
    _random_hisat_bam(bam, seed=n_parts)
    ref_stats, ref_introns = tiles.scan_batch_bam(bam, tile_size=5000)
    with ThreadPoolExecutor(4) as ex:
        stats, introns = tiles.scan_batch_bam_parallel(bam, 5000, ex, n_parts=n_parts)
    assert stats.umr_counts == ref_stats.umr_counts
    assert stats.n_reads == ref_stats.n_reads
    assert stats.n_spliced == ref_stats.n_spliced
    assert introns.counts == ref_introns.counts


@requires_pysam
def test_scan_batch_bam_parallel_in_spawned_processes(tmp_path: Path):
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor

    bam = tmp_path / "b.bam"
    _random_hisat_bam(bam, seed=42, n_reads=2000)
    ref_stats, ref_introns = tiles.scan_batch_bam(bam, tile_size=5000)
    with ProcessPoolExecutor(2, mp_context=mp.get_context("spawn")) as ex:
        stats, introns = tiles.scan_batch_bam_parallel(bam, 5000, ex, n_parts=8)
    assert stats.umr_counts == ref_stats.umr_counts
    assert (stats.n_reads, stats.n_spliced) == (ref_stats.n_reads, ref_stats.n_spliced)
    assert introns.counts == ref_introns.counts


def test_split_regions_covers_genome_once():
    refs = [("a", 23_000), ("b", 1_000), ("c", 51_000)]
    for n in (1, 2, 5, 40):
        groups = tiles.split_regions(refs, n, 5000)
        regs = [r for g in groups for r in g]
        for chrom, ln in refs:
            spans = sorted((s, e) for c, s, e in regs if c == chrom)
            assert spans[0][0] == 0 and spans[-1][1] == ln
            assert all(a[1] == b[0] for a, b in zip(spans, spans[1:]))
        assert all(s % 5000 == 0 for _, s, _ in regs)
        assert len(groups) <= n + 1

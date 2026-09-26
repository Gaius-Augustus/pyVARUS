"""Tests for varus.align.

The full hisat2 | samtools sort pipe is exercised on the cluster, not here. We
unit-test the log parser and verify that the wrapper checks for the binaries.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from varus import align


def test_parse_hisat2_log_paired(tmp_path: Path):
    # Real-shape HISAT2 paired-end summary; the legacy regex picks up the
    # "aligned concordantly exactly 1 time" line.
    log = tmp_path / "Log.final.out"
    log.write_text(
        "50000 reads; of these:\n"
        "  50000 (100.00%) were paired; of these:\n"
        "    21662 (43.32%) aligned concordantly exactly 1 time\n"
        "    100 (0.20%) aligned concordantly >1 times\n"
        "    1234 (2.47%) aligned discordantly 1 time\n"
    )
    stats = align.parse_hisat2_log(log, batch_size=50000)
    assert stats["num_uniq"] == 21662
    assert abs(stats["uniq_pct"] - 43.324) < 0.01


def test_parse_hisat2_log_single_end(tmp_path: Path):
    log = tmp_path / "Log.final.out"
    log.write_text(
        "50000 reads; of these:\n"
        "  29826 (59.65%) aligned exactly 1 time\n"
        "  900 (1.80%) aligned >1 times\n"
    )
    stats = align.parse_hisat2_log(log, batch_size=50000)
    assert stats["num_uniq"] == 29826


def test_parse_hisat2_log_no_match(tmp_path: Path):
    log = tmp_path / "Log.final.out"
    log.write_text("nothing useful here\n")
    stats = align.parse_hisat2_log(log, batch_size=50000)
    assert stats["num_uniq"] == 0
    assert stats["uniq_pct"] == 0.0


def test_parse_hisat2_log_missing_file(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        align.parse_hisat2_log(tmp_path / "nope.out", batch_size=10)


def test_align_batch_requires_hisat2(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(align.shutil, "which", lambda t: None)
    with pytest.raises(RuntimeError, match="hisat2 not found"):
        align.align_batch_hisat2(
            r1=tmp_path / "r1.fa",
            r2=None,
            index_prefix=tmp_path / "idx",
            batch_dir=tmp_path / "b",
        )


def test_align_batch_requires_samtools(tmp_path: Path, monkeypatch):
    # hisat2 found, samtools missing
    monkeypatch.setattr(align.shutil, "which",
                        lambda t: "/fake/hisat2" if t == "hisat2" else None)
    with pytest.raises(RuntimeError, match="samtools not found"):
        align.align_batch_hisat2(
            r1=tmp_path / "r1.fa",
            r2=None,
            index_prefix=tmp_path / "idx",
            batch_dir=tmp_path / "b",
        )


# ---------------------------------------------------------------------------
# minimap2 long-read path
# ---------------------------------------------------------------------------

def test_preset_for_platform_known():
    assert align.preset_for_platform("PACBIO_SMRT") == "pacbio"
    assert align.preset_for_platform("OXFORD_NANOPORE") == "ont"


def test_preset_for_platform_case_insensitive():
    assert align.preset_for_platform("pacbio_smrt") == "pacbio"
    assert align.preset_for_platform("  oxford_nanopore  ") == "ont"


def test_preset_for_platform_unknown_defaults_to_pacbio(caplog):
    import logging
    with caplog.at_level(logging.WARNING, logger="varus.align"):
        assert align.preset_for_platform("ILLUMINA") == "pacbio"
        assert align.preset_for_platform("") == "pacbio"
    assert any("Unknown" in rec.message for rec in caplog.records)


def test_align_batch_minimap2_unknown_preset(tmp_path: Path):
    with pytest.raises(ValueError, match="unknown long-read preset"):
        align.align_batch_minimap2(
            reads=tmp_path / "r.fa",
            index=tmp_path / "g.mmi",
            batch_dir=tmp_path / "b",
            preset="bogus",
        )


def test_align_batch_minimap2_requires_minimap2(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(align.shutil, "which", lambda t: None)
    with pytest.raises(RuntimeError, match="minimap2 not found"):
        align.align_batch_minimap2(
            reads=tmp_path / "r.fa",
            index=tmp_path / "g.mmi",
            batch_dir=tmp_path / "b",
        )


def test_align_batch_minimap2_requires_samtools(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(
        align.shutil, "which",
        lambda t: "/fake/minimap2" if t == "minimap2" else None,
    )
    with pytest.raises(RuntimeError, match="samtools not found"):
        align.align_batch_minimap2(
            reads=tmp_path / "r.fa",
            index=tmp_path / "g.mmi",
            batch_dir=tmp_path / "b",
        )


def test_align_batch_minimap2_command_construction(tmp_path: Path, monkeypatch):
    """Capture the cmd lists handed to subprocess.Popen."""
    monkeypatch.setattr(
        align.shutil, "which",
        lambda t: f"/fake/{t}",
    )
    captured: list[list[str]] = []

    class FakeProc:
        def __init__(self, cmd, **kwargs):
            captured.append(cmd)
            self.stdout = kwargs.get("stdin") and None
            # mm2 needs a stdout pipe for the sort to read; fake one
            if "stdout" in kwargs and kwargs["stdout"] is align.subprocess.PIPE:
                # provide a closeable object
                import io
                self.stdout = io.BytesIO()
            self._returncode = 0

        def wait(self):
            return self._returncode

    monkeypatch.setattr(align.subprocess, "Popen", FakeProc)

    junc = tmp_path / "junc.bed"
    junc.write_text("chr1\t0\t10\tj\t1\t+\t0\t10\t0\t2\t1,1\t0,9\n")

    align.align_batch_minimap2(
        reads=tmp_path / "r.fa",
        index=tmp_path / "g.mmi",
        batch_dir=tmp_path / "b",
        threads=8,
        preset="ont",
        junc_bed=junc,
    )

    assert len(captured) == 2  # mm2 then sort
    mm2_cmd, sort_cmd = captured
    assert mm2_cmd[0] == "minimap2"
    assert mm2_cmd[1:3] == ["-t", "8"]
    # ONT preset: -ax splice -uf -k14
    assert "-ax" in mm2_cmd and "splice" in mm2_cmd
    assert "-uf" in mm2_cmd and "-k14" in mm2_cmd
    assert "--junc-bed" in mm2_cmd
    assert str(junc) in mm2_cmd
    # Reads file is the last positional arg.
    assert mm2_cmd[-1] == str(tmp_path / "r.fa")
    assert mm2_cmd[-2] == str(tmp_path / "g.mmi")
    assert sort_cmd[0] == "samtools"
    assert "sort" in sort_cmd


def test_align_batch_minimap2_pacbio_preset(tmp_path: Path, monkeypatch):
    """PacBio preset is just '-ax splice', no -uf / -k14."""
    monkeypatch.setattr(align.shutil, "which", lambda t: f"/fake/{t}")
    captured: list[list[str]] = []

    class FakeProc:
        def __init__(self, cmd, **kwargs):
            captured.append(cmd)
            import io
            self.stdout = io.BytesIO() if kwargs.get("stdout") is align.subprocess.PIPE else None

        def wait(self):
            return 0

    monkeypatch.setattr(align.subprocess, "Popen", FakeProc)

    align.align_batch_minimap2(
        reads=tmp_path / "r.fa",
        index=tmp_path / "g.mmi",
        batch_dir=tmp_path / "b",
        preset="pacbio",
    )

    mm2_cmd = captured[0]
    assert "-uf" not in mm2_cmd
    assert "-k14" not in mm2_cmd
    assert "splice" in mm2_cmd


# count_minimap2_quality is gated on pysam since it scans a real BAM.
try:
    from tests.conftest import requires_pysam
except ImportError:
    requires_pysam = pytest.mark.skip(reason="conftest not found")


@requires_pysam
def test_count_minimap2_quality(tmp_path: Path):
    """Build a tiny BAM and check primary/MAPQ accounting."""
    import pysam

    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"LN": 1000, "SN": "chr1"}],
    }
    bam_path = tmp_path / "test.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        # 3 primary, MAPQ=60
        for i in range(3):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.flag = 0
            a.reference_id = 0
            a.reference_start = i * 10
            a.mapping_quality = 60
            a.cigar = [(0, 10)]
            a.query_sequence = "A" * 10
            a.query_qualities = pysam.qualitystring_to_array("I" * 10)
            bam.write(a)
        # 1 primary MAPQ=0
        a = pysam.AlignedSegment()
        a.query_name = "r_low"
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 0
        a.cigar = [(0, 10)]
        a.query_sequence = "A" * 10
        a.query_qualities = pysam.qualitystring_to_array("I" * 10)
        bam.write(a)
        # 1 secondary
        a = pysam.AlignedSegment()
        a.query_name = "r_sec"
        a.flag = 256
        a.reference_id = 0
        a.reference_start = 200
        a.mapping_quality = 60
        a.cigar = [(0, 10)]
        a.query_sequence = "A" * 10
        a.query_qualities = pysam.qualitystring_to_array("I" * 10)
        bam.write(a)
        # 1 unmapped
        a = pysam.AlignedSegment()
        a.query_name = "r_unmap"
        a.flag = 4
        a.reference_id = -1
        bam.write(a)

    pysam.index(str(bam_path))
    stats = align.count_minimap2_quality(bam_path, min_mapq=1)
    # 4 primaries (3 high-MAPQ + 1 MAPQ=0), 3 uniques.
    assert stats["num_uniq"] == 3.0
    assert abs(stats["uniq_pct"] - 75.0) < 1e-6


def test_count_minimap2_quality_missing_file(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        align.count_minimap2_quality(tmp_path / "nope.bam")


# ---------------------------------------------------------------------------
# minimap2 multi-part indexes (genomes > ~8 Gbp)
# ---------------------------------------------------------------------------


def _write_mmi(path: Path, parts, b: int = 4, flag: int = 0) -> list:
    """Write a synthetic .mmi in minimap2's mm_idx_dump layout; ``parts`` is
    a list of lists of (name, length). Returns the byte size of every part."""
    import struct
    sizes = []
    with open(path, "wb") as fh:
        for pi, seqs in enumerate(parts):
            start = fh.tell()
            fh.write(b"MMI\x02" + struct.pack("<5I", 10, 15, b, len(seqs), flag))
            for name, ln in seqs:
                fh.write(bytes([len(name)]) + name.encode() + struct.pack("<I", ln))
            for i in range(1 << b):                    # bucket i: i positions, i % 3 hash pairs
                fh.write(struct.pack("<I", i) + b"\x01" * 8 * i)
                fh.write(struct.pack("<I", i % 3) + b"\x02" * 16 * (i % 3))
            if not flag & 0x2:
                fh.write(b"\x03" * 4 * ((sum(ln for _, ln in seqs) + 7) // 8))
            sizes.append(fh.tell() - start)
    return sizes


def test_minimap2_index_parts_follow_the_dump_layout(tmp_path: Path):
    one = tmp_path / "one.mmi"
    assert align.minimap2_index_parts(one) is None                     # missing
    sizes = _write_mmi(one, [[("chr1", 1000), ("chr2", 37)]])
    assert align.minimap2_index_parts(one) == sizes
    split = tmp_path / "split.mmi"
    sizes = _write_mmi(split, [[("chr1", 1000)], [("chr2", 5000), ("", 3)], [("c3", 9)]])
    assert align.minimap2_index_parts(split) == sizes
    assert sum(sizes) == split.stat().st_size
    noseq = tmp_path / "noseq.mmi"
    assert align.minimap2_index_parts(noseq) is None
    assert align.minimap2_index_parts(tmp_path / "noseq.mmi") is None
    sizes = _write_mmi(noseq, [[("chr1", 1000)], [("chr2", 50)]], flag=0x2)
    assert align.minimap2_index_parts(noseq) == sizes
    trunc = tmp_path / "trunc.mmi"
    trunc.write_bytes(split.read_bytes()[:-5])
    assert align.minimap2_index_parts(trunc) is None
    fasta = tmp_path / "g.fa"
    fasta.write_text(">chr1\nACGT\n")
    assert align.minimap2_index_parts(fasta) is None


def test_minimap2_split_prefix_only_for_split_indexes(tmp_path: Path, caplog):
    one, split, bad = tmp_path / "one.mmi", tmp_path / "split.mmi", tmp_path / "bad.mmi"
    _write_mmi(one, [[("chr1", 1000)]])
    _write_mmi(split, [[("chr1", 1000)], [("chr2", 1000)]])
    bad.write_bytes(b"MMI\x02" + b"\x00" * 3)                # unreadable layout
    fasta = tmp_path / "g.fa"
    fasta.write_text(">chr1\nACGT\n")
    pre = tmp_path / "w" / "x.split"
    assert align.minimap2_split_prefix(one, pre) is None
    assert align.minimap2_split_prefix(fasta, pre) is None     # small FASTA: one part
    assert align.minimap2_split_prefix(tmp_path / "missing.mmi", pre) is None
    assert align.minimap2_split_prefix(split, pre) == pre and pre.parent.is_dir()
    with caplog.at_level("WARNING", logger=align.log.name):
        assert align.minimap2_split_prefix(bad, pre) == pre     # correct for any index
    assert "part layout" in caplog.text


def test_contig_alignment_passes_split_prefix_for_split_index(tmp_path: Path, monkeypatch):
    split, one = tmp_path / "split.mmi", tmp_path / "one.mmi"
    _write_mmi(split, [[("chr1", 1000)], [("chr2", 1000)]])
    _write_mmi(one, [[("chr1", 1000)]])
    monkeypatch.setattr(align.shutil, "which", lambda t: f"/fake/{t}")
    cmds = []

    class _P:
        def __init__(self, cmd, **kw):
            cmds.append(list(cmd))

    monkeypatch.setattr(align.subprocess, "Popen", _P)
    for idx in (one, split):
        align.start_contig_alignment([tmp_path / "a.fa"], index=idx, stdout=1,
                                     log_path=tmp_path / "c" / f"{idx.stem}.minimap2.err")
    assert "--split-prefix" not in cmds[0]
    sp = cmds[1][cmds[1].index("--split-prefix") + 1]
    assert sp == str(tmp_path / "c" / "split.minimap2.split")
    assert cmds[1].index("--split-prefix") < cmds[1].index(str(split))   # before the index


def test_split_queries_concatenates_only_for_split_calls(tmp_path: Path):
    a, b = tmp_path / "a.fa", tmp_path / "b.fa"
    a.write_text(">a_0\nACGT\n>a_1\nGGCC\n")
    b.write_text(">b_0\nTTTT")                                       # no final newline
    assert align.split_queries([a, b], None) == [a, b]
    assert align.split_queries([a], tmp_path / "x.split") == [a]
    (out,) = align.split_queries([a, b, a], tmp_path / "x.split")
    assert out == tmp_path / "x.split.queries.fa"
    assert out.read_text() == ">a_0\nACGT\n>a_1\nGGCC\n>b_0\nTTTT\n>a_0\nACGT\n>a_1\nGGCC\n"
    (tmp_path / "x.split.0000.tmp").write_text("")
    align.remove_split_files(tmp_path / "x.split")
    assert sorted(p.name for p in tmp_path.iterdir()) == ["a.fa", "b.fa"]


@pytest.mark.skipif(not (__import__("shutil").which("minimap2") and __import__("shutil").which("samtools")),
                    reason="minimap2/samtools not on PATH")
def test_split_index_gives_single_index_alignments(tmp_path: Path):
    """A genome indexed in 4 parts aligns like one indexed in 1 part (no ties)."""
    import os
    import random
    import subprocess
    rng = random.Random(1)
    rnd = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    # two query files with different record counts (Logan: one per run)
    genome, queries, q2 = tmp_path / "g.fa", tmp_path / "q.fa", tmp_path / "q2.fa"
    with open(genome, "w") as g, open(queries, "w") as q, open(q2, "w") as q2h:
        for c in range(4):
            exons = [rnd(300) for _ in range(3)]
            seq = rnd(50_000) + exons[0] + "GT" + rnd(1996) + "AG" + exons[1] + "GT" + \
                rnd(1996) + "AG" + exons[2] + rnd(250_000)
            g.write(f">chr{c}\n{seq}\n")
            q.write(f">tx{c}_0\n{''.join(exons)}\n")
            if c == 2:
                q2h.write(f">other_0\n{''.join(exons)}")      # no final newline
    recs = {}
    for name, extra in (("one", []), ("split", ["-I", "200k"])):
        mmi = tmp_path / f"{name}.mmi"
        subprocess.run(["minimap2", "-x", "splice", *extra, "-d", str(mmi), str(genome)],
                       check=True, capture_output=True)
        r, w = os.pipe()
        proc = align.start_contig_alignment([queries, q2], index=mmi, stdout=w, threads=2,
                                            log_path=tmp_path / name / "c.minimap2.err")
        os.close(w)
        with os.fdopen(r) as fh:
            sam = fh.read()
        assert proc.wait() == 0
        recs[name] = (sam.count("@SQ"),
                      sorted(l.split("\t")[:6] for l in sam.splitlines() if not l.startswith("@")))
        align.remove_split_files(tmp_path / name / "c.minimap2.split")
        assert not list((tmp_path / name).glob("*.split*"))
    assert len(align.minimap2_index_parts(tmp_path / "split.mmi")) == 4
    assert recs["split"] == recs["one"] and recs["one"][0] == 4
    assert len(recs["one"][1]) == 5                   # every record of both files

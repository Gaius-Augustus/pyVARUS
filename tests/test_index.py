"""Tests for varus.index.

We don't actually run hisat2-build; we mock subprocess.run and shutil.which to
verify the command construction and error paths.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from varus import index


def test_build_hisat2_index_missing_genome(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        index.build_hisat2_index(
            genome=tmp_path / "nope.fa",
            outdir=tmp_path / "idx",
        )


def test_build_hisat2_index_missing_tool(tmp_path: Path, monkeypatch):
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nACGT\n")
    monkeypatch.setattr(index.shutil, "which", lambda _: None)
    with pytest.raises(RuntimeError, match="hisat2-build not found"):
        index.build_hisat2_index(genome=genome, outdir=tmp_path / "idx")


def test_build_hisat2_index_invokes_correct_command(tmp_path: Path, monkeypatch):
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nACGT\n")
    outdir = tmp_path / "idx"

    monkeypatch.setattr(index.shutil, "which", lambda _: "/fake/hisat2-build")
    captured = {}

    def fake_run(cmd, check):
        captured["cmd"] = cmd
        captured["check"] = check
        for i in range(1, 9):
            Path(f"{cmd[-1]}.{i}.ht2").write_text("x")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(index.subprocess, "run", fake_run)

    out = index.build_hisat2_index(
        genome=genome, outdir=outdir, threads=8, prefix="myidx"
    )
    assert out == outdir / "myidx"
    assert captured["check"] is True
    assert captured["cmd"] == [
        "hisat2-build", "-p", "8", str(genome), str(outdir / "myidx.partial"),
    ]
    # Built under the temporary prefix, then renamed to the final one.
    assert sorted(f.name for f in outdir.iterdir()) == [
        f"myidx.{i}.ht2" for i in range(1, 9)]


def test_build_hisat2_index_killed_build_leaves_old_index(tmp_path: Path, monkeypatch):
    """A build killed after writing part of its files (job time limit) must
    not replace or pose as the index: hisat2-align segfaults on such files."""
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nACGT\n")
    outdir = tmp_path / "idx"
    outdir.mkdir()
    (outdir / "hisatidx.1.ht2").write_text("old")
    (outdir / "hisatidx.partial.5.ht2").write_text("stale")
    monkeypatch.setattr(index.shutil, "which", lambda _: "/fake/hisat2-build")

    def killed_run(cmd, check):
        assert not Path(f"{cmd[-1]}.5.ht2").exists()   # stale partial removed first
        Path(f"{cmd[-1]}.1.ht2").write_text("new")
        raise subprocess.CalledProcessError(-9, cmd)

    monkeypatch.setattr(index.subprocess, "run", killed_run)
    with pytest.raises(subprocess.CalledProcessError):
        index.build_hisat2_index(genome=genome, outdir=outdir)
    assert (outdir / "hisatidx.1.ht2").read_text() == "old"


def test_build_hisat2_index_replaces_large_index(tmp_path: Path, monkeypatch):
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nACGT\n")
    outdir = tmp_path / "idx"
    outdir.mkdir()
    (outdir / "hisatidx.1.ht2l").write_text("old")
    (outdir / "hisatidx2.1.ht2").write_text("other index")
    monkeypatch.setattr(index.shutil, "which", lambda _: "/fake/hisat2-build")
    monkeypatch.setattr(
        index.subprocess, "run",
        lambda cmd, check: Path(f"{cmd[-1]}.1.ht2").write_text("new")
        and subprocess.CompletedProcess(cmd, 0))
    index.build_hisat2_index(genome=genome, outdir=outdir)
    assert sorted(f.name for f in outdir.iterdir()) == ["hisatidx.1.ht2", "hisatidx2.1.ht2"]


def test_build_minimap2_index_missing_genome(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        index.build_minimap2_index(
            genome=tmp_path / "nope.fa",
            outdir=tmp_path / "idx",
        )


def test_build_minimap2_index_missing_tool(tmp_path: Path, monkeypatch):
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nACGT\n")
    monkeypatch.setattr(index.shutil, "which", lambda _: None)
    with pytest.raises(RuntimeError, match="minimap2 not found"):
        index.build_minimap2_index(genome=genome, outdir=tmp_path / "idx")


def test_build_minimap2_index_invokes_correct_command(tmp_path: Path, monkeypatch):
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nACGT\n")
    outdir = tmp_path / "idx"

    monkeypatch.setattr(index.shutil, "which", lambda _: "/fake/minimap2")
    captured = {}

    def fake_run(cmd, check):
        captured["cmd"] = cmd
        captured["check"] = check
        Path(cmd[cmd.index("-d") + 1]).write_text("x")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(index.subprocess, "run", fake_run)

    out = index.build_minimap2_index(
        genome=genome, outdir=outdir, threads=8, prefix="myidx"
    )
    assert out == outdir / "myidx.mmi"
    assert captured["check"] is True
    assert captured["cmd"] == [
        "minimap2", "-t", "8", "-x", "splice",
        "-d", str(outdir / "myidx.mmi.partial"), str(genome),
    ]
    assert [f.name for f in outdir.iterdir()] == ["myidx.mmi"]


def test_build_minimap2_index_part_bases(tmp_path: Path, monkeypatch):
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nACGT\n")
    monkeypatch.setattr(index.shutil, "which", lambda _: "/fake/minimap2")
    captured = {}
    def fake_run(cmd, check):
        captured["cmd"] = cmd
        Path(cmd[cmd.index("-d") + 1]).write_text("x")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(index.subprocess, "run", fake_run)
    index.build_minimap2_index(genome=genome, outdir=tmp_path / "idx", threads=8,
                               part_bases=14_760_000_000)
    cmd = captured["cmd"]
    assert cmd[cmd.index("-I") + 1] == "14760000000" and cmd.index("-I") < cmd.index("-d")


def test_check_sequence_lengths(tmp_path: Path, monkeypatch):
    pytest.importorskip("pyfaidx")
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\n" + "ACGT" * 30 + "\n>chr2\nACGT\n")
    index.check_sequence_lengths(genome)                 # builds g.fa.fai
    assert (tmp_path / "g.fa.fai").is_file()
    monkeypatch.setattr(index, "MAX_SEQUENCE_LENGTH", 100)
    with pytest.raises(ValueError, match=r"1 sequence\(s\) longer.*chr1"):
        index.check_sequence_lengths(genome)
    bad = tmp_path / "bad.fa.gz"
    bad.write_bytes(b"not gzip")
    index.check_sequence_lengths(bad)                    # cannot index: warns only

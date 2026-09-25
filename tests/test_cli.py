"""Smoke tests for the argparse interface."""

from __future__ import annotations

import pytest

from varus import cli


def test_parser_builds():
    p = cli.build_parser()
    assert p.prog == "varus"


def test_runlist_parser_accepts_minimal_args():
    args = cli.build_parser().parse_args(
        ["runlist", "Foo bar"]
    )
    assert args.cmd == "runlist"
    assert args.species == "Foo bar"
    assert args.paired_only is False
    assert args.max_runs == 0


def test_index_parser_accepts_minimal_args():
    args = cli.build_parser().parse_args(
        ["index", "genome.fa"]
    )
    assert args.cmd == "index"
    assert str(args.genome) == "genome.fa"
    assert args.threads == 4


def test_run_subcommand_parser_defaults(tmp_path):
    """'varus run' parser is wired up and applies correct defaults."""
    args = cli.build_parser().parse_args(
        [
            "run", "Foo bar", "genome.fa",
            "--runlist", str(tmp_path / "Runlist.tsv"),
            "--index", str(tmp_path / "genome/"),
        ]
    )
    assert args.cmd == "run"
    # Batch size and min-mapq default to None so main() can pick a mode-aware value.
    assert args.batch_size is None
    assert args.min_mapq is None
    assert args.tile_size == 5_000
    assert args.max_batches == 1_000
    assert args.threads == 4
    assert args.keep_batches is False
    assert args.bootstrap_all is False
    assert args.longreads is False


def test_index_parser_longreads_flag(tmp_path):
    """--longreads on the index subcommand toggles minimap2 mode."""
    args = cli.build_parser().parse_args(
        ["index", "genome.fa", "--longreads"]
    )
    assert args.longreads is True
    # Prefix default is None so dispatcher can pick 'mm2idx' / 'hisatidx'.
    assert args.prefix is None


def test_runlist_parser_longreads_flag():
    args = cli.build_parser().parse_args(
        ["runlist", "Foo bar", "--longreads"]
    )
    assert args.longreads is True


def test_run_parser_longreads(tmp_path):
    args = cli.build_parser().parse_args(
        [
            "run", "Foo bar", "genome.fa",
            "--runlist", str(tmp_path / "Runlist.tsv"),
            "--index", str(tmp_path / "genome/mm2idx.mmi"),
            "--longreads",
        ]
    )
    assert args.longreads is True


def test_subcommand_required():
    with pytest.raises(SystemExit):
        cli.build_parser().parse_args([])


def test_run_parser_logan_default(tmp_path):
    """Logan is on by default: no --logan-dir, --no-logan off."""
    args = cli.build_parser().parse_args(
        ["run", "Foo bar", "genome.fa",
         "--runlist", str(tmp_path / "Runlist.tsv"),
         "--index", str(tmp_path / "genome/")]
    )
    assert args.no_logan is False
    assert args.logan_dir is None
    args = cli.build_parser().parse_args(
        ["run", "Foo bar", "genome.fa",
         "--runlist", str(tmp_path / "Runlist.tsv"),
         "--index", str(tmp_path / "genome/"), "--no-logan"]
    )
    assert args.no_logan is True


def _prescreen_args(tmp_path, longreads=False):
    from types import SimpleNamespace
    args = SimpleNamespace(runlist=tmp_path / "Runlist.tsv", longreads=longreads)
    cfg = SimpleNamespace(
        genome=tmp_path / "genome.fa", outdir=tmp_path, index_prefix=tmp_path / "idx",
        threads=4, tile_size=5000, batch_size=50_000, logan_prior_batches=1.0, seed=None,
    )
    return args, cfg


def test_logan_prescreen_reuses_existing_dir(tmp_path, monkeypatch):
    """An existing <outdir>/logan/LoganRanking.tsv is used as is; nothing runs."""
    args, cfg = _prescreen_args(tmp_path)
    (tmp_path / "logan").mkdir()
    (tmp_path / "logan" / "LoganRanking.tsv").write_text("acc\tstatus\n")
    monkeypatch.setattr(cli.shutil, "which", lambda name: None)  # would fail if reached
    assert cli.logan_prescreen(args, cfg) == tmp_path / "logan"


def test_logan_prescreen_needs_minimap2(tmp_path, monkeypatch):
    args, cfg = _prescreen_args(tmp_path)
    monkeypatch.setattr(cli.shutil, "which", lambda name: None)
    with pytest.raises(SystemExit, match="--no-logan"):
        cli.logan_prescreen(args, cfg)


@pytest.mark.parametrize("rc,expect_dir", [(0, True), (3, False), (4, False)])
def test_logan_prescreen_exit_codes(tmp_path, monkeypatch, rc, expect_dir):
    """Exit 0 returns the Logan dir; 3 without a ranking and 4 fall back to no prior."""
    import varus.logan as logan_mod
    args, cfg = _prescreen_args(tmp_path, longreads=True)
    monkeypatch.setattr(cli.shutil, "which", lambda name: "/usr/bin/minimap2")
    seen = {}

    def fake_run_logan(lcfg):
        seen["cfg"] = lcfg
        return rc

    monkeypatch.setattr(logan_mod, "run_logan", fake_run_logan)
    out = cli.logan_prescreen(args, cfg)
    assert (out == tmp_path / "logan") is expect_dir
    lcfg = seen["cfg"]
    assert lcfg.outdir == tmp_path and lcfg.threads == 4 and lcfg.seed == 1
    assert lcfg.longreads is True and lcfg.mmi == tmp_path / "idx"


def test_logan_prescreen_exit3_keeps_filter(tmp_path, monkeypatch):
    """Exit 3 with a written ranking returns the dir so rejected runs are dropped."""
    import varus.logan as logan_mod
    args, cfg = _prescreen_args(tmp_path)
    monkeypatch.setattr(cli.shutil, "which", lambda name: "/usr/bin/minimap2")

    def fake_run_logan(lcfg):
        lcfg.logan_dir.mkdir(parents=True)
        (lcfg.logan_dir / "LoganRanking.tsv").write_text("acc\tstatus\nSRR1\trejected\n")
        return 3

    monkeypatch.setattr(logan_mod, "run_logan", fake_run_logan)
    assert cli.logan_prescreen(args, cfg) == tmp_path / "logan"

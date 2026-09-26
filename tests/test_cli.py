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
    assert lcfg.mmi == tmp_path / "idx"


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


def _help(argv, capsys):
    with pytest.raises(SystemExit) as e:
        cli.build_parser().parse_args(argv)
    assert e.value.code == 0
    return capsys.readouterr().out


@pytest.mark.parametrize("cmd,basic,expert", [
    ("run", ["--max-batches", "--no-logan", "--longreads"],
     ["--parallel-downloads", "--tile-size", "--logan-keep-unprocessed", "--advanced"]),
    ("logan", ["--max-candidates", "--mmi"],
     ["--min-tiles-frac", "--align-groups", "--select-top"]),
])
def test_help_hides_expert_options(capsys, cmd, basic, expert):
    """--help lists the everyday options; --help-all adds the expert groups."""
    short = _help([cmd, "--help"], capsys)
    full = _help([cmd, "--help-all"], capsys)
    for opt in basic:
        assert opt in short and opt in full
    for opt in expert:
        assert opt not in short and opt in full
    assert "expert:" not in short and "expert:" in full
    assert "--help-all" in short


def test_expert_options_still_parse(tmp_path):
    args = cli.build_parser().parse_args(
        ["run", "Foo bar", "genome.fa", "--runlist", str(tmp_path / "r.tsv"),
         "--index", str(tmp_path / "idx"), "--parallel-downloads", "3",
         "--logan-keep-unprocessed", "--advanced", "lambda=1"]
    )
    assert args.parallel_downloads == 3 and args.logan_keep_unprocessed
    assert args.advanced == ["lambda=1"]


@pytest.mark.parametrize("flag", ["--prefetch", "--no-align-ahead", "--no-hisat2-mm", "--logan-only",
                                  "--keep-unaligned", "--logan-prior-first-only",
                                  "--no-logan-seed-db", "--no-logan-bootstrap"])
def test_removed_run_options_are_rejected(tmp_path, flag):
    with pytest.raises(SystemExit):
        cli.build_parser().parse_args(
            ["run", "Foo bar", "genome.fa", "--runlist", str(tmp_path / "r.tsv"),
             "--index", str(tmp_path / "idx"), flag]
        )


def _run_main_with_logan(tmp_path, monkeypatch, argv_extra, statuses):
    """Drive `varus run` through cli.main with the loop and the pre-screen
    mocked; return the runs handed to the Controller."""
    import varus.controller as ctl
    import varus.logan as logan_mod
    from types import SimpleNamespace
    from varus.runlist import RunRecord

    rl = tmp_path / "Runlist.tsv"
    rl.write_text("")
    ldir = tmp_path / "logan"
    ldir.mkdir()
    (ldir / "LoganRanking.tsv").write_text("acc\tstatus\n")  # -> pre-screen is reused

    def fake_load_runs(path, batch_size, rng):
        recs = [RunRecord(accession=a, total_spots=n, total_bases=n * 100, avg_len=100.0,
                          paired=False, colorspace=False, platform="ILLUMINA", bioproject="")
                for a, n in (("A", 1_000_000), ("B", 500_000), ("U", 5_000_000))]
        return [ctl.RunState.from_record(r, batch_size, rng) for r in recs]

    seen = {}

    class FakeController:
        def __init__(self, cfg, runs, logan=None):
            seen["runs"] = [r.record.accession for r in runs]
            seen["cfg"] = cfg
        def run(self):
            return 0

    monkeypatch.setattr(ctl, "load_runs", fake_load_runs)
    monkeypatch.setattr(ctl, "Controller", FakeController)
    monkeypatch.setattr(logan_mod, "load_logan", lambda d: SimpleNamespace(
        status=statuses, rank={"A": 1, "B": 2}, tiles={}, yield_pct={}, introns=None,
        splice_sites=None, junc_bed=None, params={}, counts={}, acceptance_rate=None))
    rc = cli.main(["run", "Foo bar", str(tmp_path / "g.fa"), "--runlist", str(rl),
                   "--index", str(tmp_path / "idx"), "--outdir", str(tmp_path)] + argv_extra)
    assert rc == 0
    return seen


def test_run_drops_unprocessed_runs_by_default(tmp_path, monkeypatch):
    """A and B hold 30 batches: enough for --max-batches 30, so U is dropped."""
    seen = _run_main_with_logan(tmp_path, monkeypatch, ["--max-batches", "30"],
                                {"A": "accepted", "B": "accepted"})
    assert seen["runs"] == ["A", "B"]


def test_run_expands_to_unprocessed_runs_when_capacity_is_short(tmp_path, monkeypatch):
    """With --max-batches 1000 the accepted runs cannot fill the run: U stays."""
    seen = _run_main_with_logan(tmp_path, monkeypatch, [],
                                {"A": "accepted", "B": "accepted"})
    assert seen["runs"] == ["A", "B", "U"]
    assert seen["cfg"].max_batches == 1000


def test_run_keep_unprocessed_flag(tmp_path, monkeypatch):
    seen = _run_main_with_logan(tmp_path, monkeypatch,
                                ["--max-batches", "30", "--logan-keep-unprocessed"],
                                {"A": "accepted", "B": "accepted"})
    assert seen["runs"] == ["A", "B", "U"]


def test_run_keeps_unprocessed_when_nothing_accepted(tmp_path, monkeypatch):
    seen = _run_main_with_logan(tmp_path, monkeypatch, ["--max-batches", "5"],
                                {"A": "rejected"})
    assert seen["runs"] == ["B", "U"]

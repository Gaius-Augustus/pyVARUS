#!/usr/bin/env python3
"""Benchmark pyVARUS speed-ups against a baseline run.

Two modes:

``prepare``
    Write one SLURM job script per benchmark arm for a genome + runlist
    (used on the cluster; each job runs ``varus [logan] run`` with the arm's
    flags inside the pyVARUS container or a conda env and records wall time).

``report``
    Parse the ``BatchTimings.tsv`` / ``Coverage.csv`` / ``introns.gff`` of
    every finished arm (plus an optional legacy log for the baseline) and
    print a Markdown table: wall time per phase, rejected batches, final
    score S = Σ log1p(c_j), tiles covered, intron count and Jaccard vs the
    baseline, time-to-90 %-of-baseline-S.

Arms (see ``ARMS``): baseline flags, in-loop speed-ups, +parallel downloads,
+prefetch, +Logan prior (1000 and 500 batches), +profit condition. Everything
is driven by the same seed so pick sequences are comparable where the
algorithm is unchanged.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Arms
# ---------------------------------------------------------------------------

@dataclass
class Arm:
    name: str
    run_flags: str
    logan: bool = False
    max_batches: int = 1000
    note: str = ""


ARMS: List[Arm] = [
    Arm("A0_baseline_flags", "--merge-every 0 --no-hisat2-mm --keep-unaligned",
        note="new code, legacy behaviour (serial, single final merge)"),
    Arm("A1_inloop", "", note="incremental DB, --mm, --no-unal, rolling merge"),
    Arm("A2_parallel3", "--parallel-downloads 3"),
    Arm("A3_parallel3_prefetch", "--parallel-downloads 3 --prefetch"),
    Arm("A4_logan_1000", "--parallel-downloads 3 --prefetch", logan=True),
    Arm("A5_logan_500", "--parallel-downloads 3 --prefetch", logan=True, max_batches=500),
    Arm("A6_logan_profit", "--parallel-downloads 3 --prefetch --profit-condition", logan=True),
]

SLURM_TEMPLATE = """#!/bin/bash -l
#SBATCH --job-name=varus_{arm}_{tag}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --partition={partition}
#SBATCH --output={outdir}/{arm}.slurm.out
# login shell (-l) so `module` is available inside the job
set -euo pipefail
{env}
OUT={outdir}/{arm}
mkdir -p "$OUT"
cd "$OUT"
cp -n {runlist} Runlist.tsv
S0=$(date +%s)
{logan_cmd}
S1=$(date +%s)
/usr/bin/time -p -o runtime.varus.txt \\
  varus run "{species}" {genome} \\
    --runlist {run_runlist} --index {index} --outdir . \\
    --threads {cpus} --seed {seed} --max-batches {max_batches} {flags} {logan_flags} \\
    > varus.log 2>&1 || echo "varus exit $?" >> varus.log
S2=$(date +%s)
echo "logan_seconds=$((S1-S0))" > phases.txt
echo "run_seconds=$((S2-S1))" >> phases.txt
"""


def cmd_prepare(a: argparse.Namespace) -> int:
    outdir = Path(a.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    env = a.env_setup or ""
    arms = [x for x in ARMS if not a.arms or x.name in a.arms]
    for arm in arms:
        logan_cmd = ""
        logan_flags = ""
        run_runlist = "Runlist.tsv"
        if arm.logan:
            logan_cmd = (
                f'varus logan {a.genome} --runlist Runlist.tsv --outdir . '
                f'--threads {a.cpus} --max-candidates {a.logan_max_candidates} '
                f'--select-top {a.logan_select_top} > logan.log 2>&1 || echo "logan exit $?" >> logan.log'
            )
            logan_flags = "--logan-dir logan"
            run_runlist = "Runlist.logan.tsv"
        script = SLURM_TEMPLATE.format(
            arm=arm.name, tag=a.tag, cpus=a.cpus, mem=a.mem, time=a.time,
            partition=a.partition, outdir=outdir, env=env, runlist=Path(a.runlist).resolve(),
            species=a.species, genome=Path(a.genome).resolve(), index=Path(a.index).resolve(),
            seed=a.seed, max_batches=arm.max_batches, flags=arm.run_flags,
            logan_cmd=logan_cmd, logan_flags=logan_flags, run_runlist=run_runlist,
        )
        path = outdir / f"{arm.name}.sbatch"
        path.write_text(script)
        print(path)
    return 0


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

@dataclass
class ArmResult:
    name: str
    n_batches: int = 0
    n_rejected: int = 0
    wall_loop: float = 0.0
    t_download: float = 0.0
    t_align: float = 0.0
    t_scan: float = 0.0
    t_db: float = 0.0
    t_est: float = 0.0
    logan_seconds: float = 0.0
    run_seconds: float = 0.0
    score: float = 0.0
    tiles1: int = 0
    tiles10: int = 0
    n_introns: int = 0
    introns: set = field(default_factory=set)
    score_curve: List[Tuple[float, float]] = field(default_factory=list)  # (cum seconds, S)
    exit_note: str = ""


def _read_timings(path: Path) -> List[dict]:
    rows = []
    with path.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append(r)
    return rows


def _read_coverage(path: Path) -> Dict[str, int]:
    cov = {}
    with path.open() as f:
        next(f, None)
        for line in f:
            if ";" not in line:
                continue
            k, v = line.rstrip("\n").split(";")
            cov[k] = int(v)
    return cov


def _read_introns(path: Path) -> set:
    s = set()
    with path.open() as f:
        for line in f:
            p = line.split("\t")
            if len(p) >= 5 and p[2] == "intron":
                s.add((p[0], p[3], p[4]))
    return s


def parse_arm(d: Path) -> ArmResult:
    res = ArmResult(name=d.name)
    t = d / "BatchTimings.tsv"
    if t.is_file():
        rows = _read_timings(t)
        res.n_batches = len(rows)
        res.n_rejected = sum(1 for r in rows if r["success"] == "0")
        for k, attr in (("t_download", "t_download"), ("t_align", "t_align"),
                        ("t_scan", "t_scan"), ("t_db", "t_db"), ("t_estimate", "t_est")):
            setattr(res, attr, sum(float(r[k]) for r in rows))
        res.wall_loop = sum(float(r["t_batch"]) for r in rows)
    ph = d / "phases.txt"
    if ph.is_file():
        for line in ph.read_text().splitlines():
            k, _, v = line.partition("=")
            if k == "logan_seconds":
                res.logan_seconds = float(v)
            elif k == "run_seconds":
                res.run_seconds = float(v)
    cov = d / "Coverage.csv"
    if cov.is_file():
        c = _read_coverage(cov)
        res.score = sum(math.log1p(v) for v in c.values())
        res.tiles1 = sum(1 for v in c.values() if v >= 1)
        res.tiles10 = sum(1 for v in c.values() if v >= 10)
    gff = d / "introns.gff"
    if gff.is_file():
        res.introns = _read_introns(gff)
        res.n_introns = len(res.introns)
    log = d / "varus.log"
    if log.is_file():
        m = re.findall(r"varus exit (\d+)", log.read_text(errors="replace"))
        if m:
            res.exit_note = f"exit {m[-1]}"
    return res


_TS = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d) \w+ [\w.]+: (.*)")


def parse_legacy_log(path: Path, name: str = "baseline_log") -> ArmResult:
    """Baseline from a production varus.log (no BatchTimings.tsv)."""
    res = ArmResult(name=name)
    ev = []
    for line in path.open(errors="replace"):
        m = _TS.match(line)
        if m:
            ev.append((datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S"), m.group(2)))
    b = [i for i, (_, m) in enumerate(ev) if m.startswith("Batch ")]
    if not b:
        return res
    res.n_batches = len(b)
    res.n_rejected = sum(1 for _, m in ev if "bad quality" in m)
    res.wall_loop = (ev[b[-1]][0] - ev[b[0]][0]).total_seconds()
    for x, y in zip(b, b[1:]):
        seg = ev[x:y]
        tal = next((t for t, m in seg if m.startswith("HISAT2")), None)
        tbam = next((t for t, m in seg if "UMRs across" in m), None)
        tdl = next((t for t, m in seg if m.startswith("fastq-dump")), None)
        if tal and tbam and tdl:
            res.t_align += (tbam - tal).total_seconds()
            res.t_scan += (tdl - tbam).total_seconds()
            res.t_download += (ev[y][0] - tdl).total_seconds()
    res.run_seconds = res.wall_loop
    return res


def _fmt_h(sec: float) -> str:
    return f"{sec / 3600:.2f} h" if sec >= 3600 else f"{sec / 60:.1f} min"


def cmd_report(a: argparse.Namespace) -> int:
    root = Path(a.outdir)
    arms: List[ArmResult] = []
    if a.baseline_log:
        arms.append(parse_legacy_log(Path(a.baseline_log)))
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        if (d / "BatchTimings.tsv").is_file() or (d / "Coverage.csv").is_file():
            arms.append(parse_arm(d))
    if not arms:
        print("no arms found", file=sys.stderr)
        return 1
    base = next((x for x in arms if x.name.startswith("baseline") or x.name.startswith("A0")), arms[0])
    lines = [
        "| arm | batches | rejected | loop wall | logan | dl | align | scan | db+est | S | tiles≥1 | tiles≥10 | introns | Jaccard vs base | speedup |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for x in arms:
        jac = ""
        if base.introns and x.introns:
            jac = f"{len(base.introns & x.introns) / len(base.introns | x.introns):.3f}"
        total = (x.run_seconds or x.wall_loop) + x.logan_seconds
        btotal = (base.run_seconds or base.wall_loop) + base.logan_seconds
        speed = f"{btotal / total:.2f}×" if total > 0 and btotal > 0 else ""
        lines.append(
            f"| {x.name} | {x.n_batches} | {x.n_rejected} | {_fmt_h(x.wall_loop)} | "
            f"{_fmt_h(x.logan_seconds) if x.logan_seconds else '-'} | {_fmt_h(x.t_download)} | "
            f"{_fmt_h(x.t_align)} | {_fmt_h(x.t_scan)} | {_fmt_h(x.t_db + x.t_est)} | "
            f"{x.score:.0f} | {x.tiles1} | {x.tiles10} | {x.n_introns} | {jac} | {speed} {x.exit_note} |"
        )
    out = "\n".join(lines)
    print(out)
    if a.json:
        Path(a.json).write_text(json.dumps([
            {k: (v if not isinstance(v, set) else len(v)) for k, v in vars(x).items()
             if k != "score_curve"} for x in arms
        ], indent=2))
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare", help="write one sbatch script per arm")
    p.add_argument("--species", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--runlist", required=True)
    p.add_argument("--index", required=True, help="HISAT2 index prefix")
    p.add_argument("--outdir", required=True)
    p.add_argument("--tag", default="bench")
    p.add_argument("--cpus", type=int, default=16)
    p.add_argument("--mem", default="60G")
    p.add_argument("--time", default="24:00:00")
    p.add_argument("--partition", default="batch")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--arms", nargs="*", default=[])
    p.add_argument("--logan-max-candidates", type=int, default=500)
    p.add_argument("--logan-select-top", type=int, default=50)
    p.add_argument("--env-setup", default="",
                   help="shell lines that put varus/hisat2/fastq-dump/minimap2 on PATH")
    p.set_defaults(func=cmd_prepare)
    r = sub.add_parser("report", help="summarise finished arms")
    r.add_argument("--outdir", required=True)
    r.add_argument("--baseline-log", default=None, help="legacy varus.log for the baseline row")
    r.add_argument("--json", default=None)
    r.set_defaults(func=cmd_report)
    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())

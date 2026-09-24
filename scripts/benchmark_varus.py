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
    threads: Optional[int] = None  # overrides --cpus for this arm


ARMS: List[Arm] = [
    Arm("A0_baseline_flags", "--merge-every 0 --no-hisat2-mm --keep-unaligned",
        note="new code, legacy behaviour (serial, single final merge)"),
    Arm("A1_inloop", "", note="incremental DB, --mm, --no-unal, rolling merge"),
    Arm("A2_parallel3", "--parallel-downloads 3"),
    Arm("A3_parallel3_prefetch", "--parallel-downloads 3 --prefetch"),
    Arm("A3t2_parallel3_prefetch", "--parallel-downloads 3 --prefetch", threads=2,
        note="A3 at 2 threads: like-for-like with the 2-thread production baseline"),
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

# Variant for clusters with a slow shared file system: genome + index go to
# node-local RAM disk, the container image, code and all working files to
# node-local disk; results are copied back to the shared outdir by an EXIT
# trap (also on failure, scancel or time-limit TERM) and both local
# directories are removed.
SLURM_LOCAL_TEMPLATE = """#!/bin/bash -l
#SBATCH --job-name=varus_{arm}_{tag}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --partition={partition}
#SBATCH --output={outdir}/{arm}.slurm.out
#SBATCH --signal=B:TERM@900
# login shell (-l) so `module` is available inside the job
set -uo pipefail
{env}
DEST={outdir}/{arm}
JOB=${{SLURM_JOB_ID:-$$}}
WORK={scratch}/varus_{arm}_$JOB
SHM={shm}/varus_{arm}_$JOB
# Disk guards: never fill node-local /tmp or /dev/shm (shared with other jobs).
NEED_TMP_GB={need_tmp_gb}; MIN_TMP_GB={min_tmp_gb}; NEED_SHM_GB={need_shm_gb}; MIN_SHM_GB={min_shm_gb}
free_gb() {{ df -BG --output=avail "$1" | tail -1 | tr -dc 0-9; }}
if [ "$(free_gb {scratch})" -lt "$NEED_TMP_GB" ] || [ "$(free_gb {shm})" -lt "$NEED_SHM_GB" ]; then
  echo "ABORT: not enough local space on $(hostname): {scratch} $(free_gb {scratch})G (need $NEED_TMP_GB), {shm} $(free_gb {shm})G (need $NEED_SHM_GB)"
  exit 2
fi
mkdir -p "$DEST" "$WORK/run" "$SHM"
finish() {{
  rc=$?
  trap - EXIT TERM INT
  [ -n "${{WATCH:-}}" ] && kill "$WATCH" 2>/dev/null
  if [ -n "${{MAIN:-}}" ]; then
    pkill -TERM -f "$WORK" 2>/dev/null; pkill -TERM -P "$MAIN" 2>/dev/null
    kill -TERM "$MAIN" 2>/dev/null; wait "$MAIN" 2>/dev/null; sleep 2
  fi
  echo "copy-back rc=$rc $(date '+%F %T')"
  # /home and /projects are slow Ceph (~0.5-7 MB/s): copy back results only, not
  # regenerable caches (minimap2 index, per-run tile caches, .sra, contigs).
  rsync -a --exclude '*.sra' --exclude 'sra/' --exclude 'logan/contigs/' \
    --exclude 'logan/genome/' --exclude '*.mmi' --exclude '*.tiles.npz' \
    --exclude 'batches/' --exclude 'merged/' {bam_exclude}"$WORK/run/" "$DEST/" \
    || echo "WARNING: copy-back failed"
  rm -rf "$WORK" "$SHM"
  exit $rc
}}
trap finish EXIT
trap 'echo "caught TERM"; exit 143' TERM INT
echo "node $(hostname) work $WORK shm $SHM"
cp {genome} "$SHM/genome.fa"
for f in {index}.*.ht2; do cp "$f" "$SHM/hisatidx.${{f#{index}.}}"; done
{stage}
cp {runlist} "$WORK/run/Runlist.tsv"
cd "$WORK/run"
main() {{
S0=$(date +%s)
{logan_cmd}
S1=$(date +%s)
/usr/bin/time -p -o runtime.varus.txt \
  varus run "{species}" "$SHM/genome.fa" \
    --runlist {run_runlist} --index "$SHM/hisatidx" --outdir . \
    --threads {cpus} --seed {seed} --max-batches {max_batches} {flags} {logan_flags} \
    > varus.log 2>&1 || echo "varus exit $?" >> varus.log
S2=$(date +%s)
echo "logan_seconds=$((S1-S0))" > phases.txt
echo "run_seconds=$((S2-S1))" >> phases.txt
[ -s VARUS.bam ] && {samtools} flagstat -@ 4 VARUS.bam > VARUS.flagstat 2>&1
}}
main &
MAIN=$!
# watchdog: stop the run (the EXIT trap copies back and cleans) if local space runs low
( while kill -0 "$MAIN" 2>/dev/null; do
    t=$(free_gb {scratch}); m=$(free_gb {shm})
    if [ "$t" -lt "$MIN_TMP_GB" ] || [ "$m" -lt "$MIN_SHM_GB" ]; then
      echo "WATCHDOG: low local space ({scratch} ${{t}}G, {shm} ${{m}}G), stopping run"
      echo "watchdog stop: {scratch} ${{t}}G {shm} ${{m}}G" >> "$WORK/run/varus.log"
      pkill -TERM -f "$WORK" ; kill -TERM "$MAIN"; break
    fi
    sleep 30
  done ) &
WATCH=$!
wait "$MAIN"
"""


def cmd_prepare(a: argparse.Namespace) -> int:
    outdir = Path(a.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    env = a.env_setup or ""
    arms = [x for x in ARMS if not a.arms or x.name in a.arms]
    stage = ""
    for item in a.stage:
        src, _, var = item.partition(":")
        src_p = Path(src).resolve()
        dst = f'"$WORK/{src_p.name}"'
        cp = "cp -r" if src_p.is_dir() else "cp"
        stage += f"{cp} {src_p} {dst}\n"
        if var:
            stage += f"export {var}={dst}\n"
    for arm in arms:
        cpus = arm.threads or a.cpus
        logan_cmd = ""
        logan_flags = ""
        run_runlist = "Runlist.tsv"
        if arm.logan:
            logan_cmd = (
                f'varus logan {a.genome} --runlist Runlist.tsv --outdir . '
                f'--threads {cpus} --max-candidates {a.logan_max_candidates} '
                f'--select-top {a.logan_select_top} > logan.log 2>&1 || echo "logan exit $?" >> logan.log'
            )
            logan_flags = "--logan-dir logan"
            run_runlist = "Runlist.logan.tsv"
        genome_arg = '"$SHM/genome.fa"' if a.local_scratch else Path(a.genome).resolve()
        if arm.logan:
            logan_cmd = logan_cmd.replace(f"varus logan {a.genome} ", f"varus logan {genome_arg} ")
        tpl = SLURM_LOCAL_TEMPLATE if a.local_scratch else SLURM_TEMPLATE
        bam_exclude = "" if arm.name in a.copy_bam_arms else "--exclude 'VARUS.bam*' "
        flags = " ".join(x for x in (arm.run_flags, a.extra_run_flags) if x)
        script = tpl.format(
            scratch=a.local_scratch, shm=a.shm, stage=stage.rstrip("\n"),
            bam_exclude=bam_exclude, samtools=a.samtools,
            need_tmp_gb=a.need_tmp_gb, min_tmp_gb=a.min_tmp_gb,
            need_shm_gb=a.need_shm_gb, min_shm_gb=a.min_shm_gb,
            arm=arm.name, tag=a.tag, cpus=cpus, mem=a.mem, time=a.time,
            partition=a.partition, outdir=outdir, env=env, runlist=Path(a.runlist).resolve(),
            species=a.species, genome=Path(a.genome).resolve(), index=Path(a.index).resolve(),
            seed=a.seed, max_batches=arm.max_batches, flags=flags,
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


def _read_intron_mult(path: Path) -> Dict[Tuple[str, str, str], int]:
    """introns.gff -> {(chrom, start, end): multiplicity}."""
    d: Dict[Tuple[str, str, str], int] = {}
    with path.open() as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 9 and p[2] == "intron":
                m = re.search(r"mult=(\d+)", p[8])
                d[(p[0], p[3], p[4])] = int(m.group(1)) if m else int(float(p[5] or 1))
    return d


INTRON_MIN_LEN, INTRON_MAX_LEN = 32, 350_000  # bam2hints defaults


def read_coding_introns(gff: Path) -> set:
    """Introns between consecutive CDS parts of each transcript in a GFF3/GTF.

    The VARUS paper's reference set ("introns in the protein-coding regions
    of genes"). CDS parts are grouped by ``Parent=`` (GFF3) or
    ``transcript_id`` (GTF).
    """
    parts: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
    with gff.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            m = re.search(r"Parent=([^;]+)", p[8]) or re.search(r'transcript_id "([^"]+)"', p[8])
            if not m:
                continue
            for parent in m.group(1).split(","):
                parts.setdefault((p[0], parent), []).append((int(p[3]), int(p[4])))
    introns = set()
    for (chrom, _), segs in parts.items():
        segs.sort()
        for (_, e1), (s2, _) in zip(segs, segs[1:]):
            if s2 - e1 > 1:
                introns.add((chrom, str(e1 + 1), str(s2 - 1)))
    return introns


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
    # production Coverage.csv / introns.gff copied next to the log as baseline_*
    cov = path.with_name("baseline_Coverage.csv")
    if cov.is_file():
        c = _read_coverage(cov)
        res.score = sum(math.log1p(v) for v in c.values())
        res.tiles1 = sum(1 for v in c.values() if v >= 1)
        res.tiles10 = sum(1 for v in c.values() if v >= 10)
    gff = path.with_name("baseline_introns.gff")
    if gff.is_file():
        res.introns = _read_introns(gff)
        res.n_introns = len(res.introns)
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
    # speedup against the production run (or A0); score and introns against A0,
    # which shares code, seed and runlist with every other arm
    base = next((x for x in arms if x.name.startswith("baseline") or x.name.startswith("A0")), arms[0])
    ref = next((x for x in arms if x.name.startswith("A0")), base)

    def jac(u: set, v: set) -> str:
        return f"{len(u & v) / len(u | v):.3f}" if u and v else ""

    lines = [
        "| arm | batches | rejected | run wall | logan | dl | align | scan | db+est | S | S/A0 "
        "| tiles≥1 | tiles≥10 | introns | Jaccard A0 | Jaccard base | speedup |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    btotal = (base.run_seconds or base.wall_loop) + base.logan_seconds
    for x in arms:
        total = (x.run_seconds or x.wall_loop) + x.logan_seconds
        speed = f"{btotal / total:.2f}×" if total > 0 and btotal > 0 else ""
        srel = f"{100 * x.score / ref.score:.1f} %" if ref.score else ""
        lines.append(
            f"| {x.name} | {x.n_batches} | {x.n_rejected} | {_fmt_h(x.run_seconds or x.wall_loop)} | "
            f"{_fmt_h(x.logan_seconds) if x.logan_seconds else '-'} | {_fmt_h(x.t_download)} | "
            f"{_fmt_h(x.t_align)} | {_fmt_h(x.t_scan)} | {_fmt_h(x.t_db + x.t_est)} | "
            f"{x.score:.0f} | {srel} | {x.tiles1} | {x.tiles10} | {x.n_introns} | "
            f"{jac(ref.introns, x.introns)} | {jac(base.introns, x.introns)} | {speed} {x.exit_note} |"
        )
    if a.ref_gff:
        ref = read_coding_introns(Path(a.ref_gff))
        # Predicted introns restricted to the bam2hints --intronsonly default
        # window (32 bp..350 kb), as in legacy VARUS and Stanke et al. 2019.
        lines += ["", f"Coding introns in {Path(a.ref_gff).name}: {len(ref)}; "
                  f"predicted introns {INTRON_MIN_LEN} bp–{INTRON_MAX_LEN // 1000} kb", "",
                  "| arm | introns | Sn | Sp | introns ≥ 2 | Sn ≥ 2 | Sp ≥ 2 "
                  "| introns ≥ 5 | Sn ≥ 5 | Sp ≥ 5 |",
                  "|---|---|---|---|---|---|---|---|---|---|"]
        for x in arms:
            src = (Path(a.baseline_log).with_name("baseline_introns.gff")
                   if x is arms[0] and a.baseline_log else root / x.name / "introns.gff")
            if not src.is_file() or not ref:
                continue
            mult = {k: v for k, v in _read_intron_mult(src).items()
                    if INTRON_MIN_LEN <= int(k[2]) - int(k[1]) + 1 <= INTRON_MAX_LEN}
            cells = [x.name]
            for t in (1, 2, 5):
                pred = {k for k, v in mult.items() if v >= t}
                tp = len(pred & ref)
                cells += [str(len(pred)), f"{tp / len(ref):.3f}",
                          f"{tp / len(pred):.3f}" if pred else ""]
            lines.append("| " + " | ".join(cells) + " |")
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
    p.add_argument("--extra-run-flags", default="",
                   help="flags appended to every arm's varus run, e.g. '--prefetch-disk-gb 40'")
    p.add_argument("--local-scratch", default="",
                   help="node-local directory (e.g. /tmp): work there, copy results back on exit")
    p.add_argument("--shm", default="/dev/shm",
                   help="RAM disk for genome + HISAT2 index in --local-scratch mode")
    p.add_argument("--need-tmp-gb", type=int, default=80,
                   help="--local-scratch: refuse to start below this much free space")
    p.add_argument("--min-tmp-gb", type=int, default=30,
                   help="--local-scratch: watchdog stops the run below this much free space")
    p.add_argument("--need-shm-gb", type=int, default=4)
    p.add_argument("--min-shm-gb", type=int, default=2)
    p.add_argument("--copy-bam-arms", nargs="*", default=[],
                   help="local-scratch mode: arms whose VARUS.bam is copied back "
                        "(default none; a flagstat is always kept)")
    p.add_argument("--samtools", default="samtools",
                   help="samtools command for the on-node flagstat, e.g. "
                        "'singularity exec $VARUS_SIF samtools'")
    p.add_argument("--stage", action="append", default=[], metavar="SRC[:ENVVAR]",
                   help="copy SRC (file or dir) to local scratch before running and export "
                        "ENVVAR=<local copy>; repeatable (e.g. the .sif and the pyVARUS checkout)")
    p.set_defaults(func=cmd_prepare)
    r = sub.add_parser("report", help="summarise finished arms")
    r.add_argument("--outdir", required=True)
    r.add_argument("--baseline-log", default=None, help="legacy varus.log for the baseline row")
    r.add_argument("--json", default=None)
    r.add_argument("--ref-gff", default=None,
                   help="reference annotation (GFF3/GTF); adds coding-intron sensitivity/specificity")
    r.set_defaults(func=cmd_report)
    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())

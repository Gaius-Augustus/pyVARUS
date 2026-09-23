# Benchmark: pyVARUS v2 speed-ups

Script: `scripts/benchmark_varus.py` (`prepare` writes one SLURM job per arm,
`report` summarises finished arms). Baseline numbers come from the two
production logs analysed on 2026-09-23 (BRAKER4 wrapper, 2 threads).

## Arms

| Arm | Flags | Question |
|---|---|---|
| A0_baseline_flags | `--merge-every 0 --no-hisat2-mm --keep-unaligned` | new code, legacy behaviour (serial, single final merge) |
| A1_inloop | (defaults) | incremental splice DB, `--mm --no-unal`, rolling merge |
| A2_parallel3 | `--parallel-downloads 3` | download overlap + in-flight-aware picks |
| A3_parallel3_prefetch | `--parallel-downloads 3 --prefetch` | local range dumps |
| A3t2_parallel3_prefetch | A3 with `--threads 2` | like-for-like against the 2-thread production baseline |
| A4_logan_1000 | A3 + `varus logan` + `--logan-dir` | prior + gate at equal batch budget |
| A5_logan_500 | A4 with `--max-batches 500` | same score with half the batches? |
| A6_logan_profit | A4 + `--profit-condition` | does the informed prior make early stopping usable? |

All arms use the same genome, `Runlist.tsv` and `--seed`. On brain the arms
run on the `snowball` partition (72 logical CPUs, 189 GB per node) with
48 CPUs and 100 GB per job, the allocation BRAKER4's `run_varus` rule
requests (`slurm_args.cpus_per_task`); only A3t2 uses 2 CPUs.

## Datasets

* Coelastrella tenuitheca GCA_051903525.1 — 9 runs, no rejections (tests
  the in-loop speed-ups in isolation; Logan can prune nothing here).
* Chlorella sorokiniana GCA_025917655.1 — 391 runs, 35 % of batches rejected
  in production (tests the Logan gate).

## Metrics (from `BatchTimings.tsv`, `Coverage.csv`, `introns.gff`, `logan/`)

* wall time per phase and total, Logan stage time;
* batches rejected;
* final score S = Σ log1p(c_j), tiles with ≥ 1 and ≥ 10 UMRs;
* intron count and Jaccard vs the A0 intron set;
* reads in `VARUS.bam`;
* Logan: candidates checked / available / accepted / rejected, coverage curve.

## Decision rules

* Adopt Logan + a lower `--max-batches` in the BRAKER4 wrapper only if A5
  reaches ≥ 95 % of A0's score and ≥ 90 % intron Jaccard in ≤ 50 % of the
  wall time.
* Adopt `--splice-db-min-mult 2` only if intron Jaccard ≥ 0.98.

## Results

_To be filled in from `scripts/benchmark_varus.py report`._

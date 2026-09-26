# Benchmark: pyVARUS v2 speed-ups

Script: `scripts/benchmark_varus.py` (`prepare` writes one SLURM job per arm,
`report` summarises finished arms). Baseline numbers come from the two
production logs analysed on 2026-09-23 (BRAKER4 wrapper, 2 threads).
[`docs/figures/speedups.svg`](figures/speedups.svg) shows how the speed-ups
fit together; the README section "Threads and machine size" lists how
`--threads` is split.

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
| A7_defaults | (v2 defaults, 2026-09-24: `--parallel-downloads 6`) | the shipped configuration without Logan |
| A8_logan_defaults | A7 + `varus logan` (16 download connections, `-K 20M`, 25-run chunks, streamed scan) + `--logan-dir` | the shipped configuration with Logan, before the sparse estimator |
| A9_defaults_sparse | A7 + sparse estimator | the shipped configuration without Logan |
| A10_logan_sparse | A8 + sparse estimator; Logan back to 8 connections, no `-K` | the shipped configuration with Logan |
| A11_logan_only | A10 + `--logan-only` | drop runs Logan never screened (A8: 183 batches, 12 of 18 rejections)? |
| A12_logan_prior | A10 + `varus logan --tile-weight ka --tile-ka-cap 0`, `varus run --logan-prior-batches 40 --logan-prior-first-only` | a prior that ranks unsampled runs (see the offline test below) |
| A9_defaults_sparse_s2/_s3 | A9 with `--seed 2` / `3` | seed-to-seed spread of S |
| A13_lam1_a01_s1–s3 | A9 + `--advanced lambda=1 pseudo-count=0.1`, seeds 1–3 | does the loop exploit super runs after one batch? |
| A14_lam3_a01_s1–s3 | A9 + `--advanced lambda=3 pseudo-count=0.1`, seeds 1–3 | same, milder smoothing |
| A15_logan_lam3_s1–s3 | A11 + `--advanced lambda=3 pseudo-count=0.1`, seeds 1–3 | Logan finds good runs early, the smoothing keeps the loop on them? |
| A16_logan_prior_lam3 | A15 + A12's uncapped tile weights and prior ×40 until the first batch | does the strong prior work once sampled runs are scored on their own counts? |
| A17_ahead_s1–s3 | A14 with align-ahead (default since 2026-09-25), seeds 1–3 | how much does aligning the next batch during the scan save? |
| A18_noahead_s1–s3 | A17 + `--no-align-ahead`, submitted together with A17 | paired control on the same code and network |
| A20_logan_merge10_s1–s3 | A22 + `--merge-batches 10` | do merged downloads remove A15's download bound? |
| A21_merge10_s1–s3 | A19 + `--merge-batches 10` | merged downloads without Logan |
| A22_logan_lam3 | A15 on the current code (align-ahead, thread budget) | control for A20 |
| A19_budget_s1–s3 | A17 with the thread budget (aligner −2 cores, −4 more during a rolling merge; `samtools sort -@ 4`) | does giving each concurrent stage its own cores cost or gain time? |
| A23_logan_pscan_s1–s3 | A20 + `--scan-workers 4` (merged batches scanned by region in 4 processes) | does the parallel scan shorten the loop, and is it lossless? |
| A24_pscan_s1–s3 | A21 + `--scan-workers 4` | the same without Logan |
| A25_logan_groups_s1–s3 | A23 + `varus logan --align-groups 3` (3 minimap2 processes per chunk, each piped into its own scanner) | is the Logan stage faster with identical rankings? |
| B1_logan_lam3_s1–s3 | current defaults + Logan (unscreened runs kept) + λ = 3, a = 0.1 | final settings, few-run species |
| B2_lam3_s1–s3 | current defaults + λ = 3, a = 0.1, no Logan | what Logan adds under λ = 3 |
| B3_logan_lam10_s1–s3 | B1 with λ = 10, a = 1 | λ = 3 against λ = 10 with Logan |
| B4_lam10_s1–s3 | current defaults, λ = 10, a = 1, no Logan | λ = 3 against λ = 10 without Logan |
| B5_loganonly_lam3_s1–s3 | B1 with `--logan-only` (many-run species) | final settings, mouse |
| B6_loganonly_lam10_s1 | B5 with λ = 10, a = 1 | λ on mouse |

"Current defaults" for B1–B6: `--parallel-downloads 6 --merge-batches 10
--scan-workers 4`, align-ahead, `varus logan --align-groups 3`, 1000
batches, 48 threads. Since λ = 3 became the default, arms without
`--advanced` get `--advanced lambda=10 pseudo-count=1` from the script.

A3/A3t2 used `--prefetch`, which was retired on 2026-09-24 and removed from
the code on 2026-09-26 together with the other experiment-only switches
(`--no-hisat2-mm`, `--keep-unaligned` in A0; `--tile-weight`,
`--tile-ka-cap`, `--logan-prior-first-only` in A12/A16; `--no-align-ahead`
in A18); those arms are kept for the record and were removed from the
script. A1/A2/A4–A6 now pass
their download count explicitly because the default became 6.

All arms use the same genome, `Runlist.tsv` and `--seed`. On brain the arms
run on the `snowball` partition (72 logical CPUs, 189 GB per node) with
48 CPUs and 100 GB per job, the allocation BRAKER4's `run_varus` rule
requests (`slurm_args.cpus_per_task`); only A3t2 uses 2 CPUs.

brain's shared `/home` file system was slow during the benchmark (a
`hisat2-build` that takes 32 s on a workstation ran for more than 20 min), so
the jobs are generated with `prepare --local-scratch /tmp --shm /dev/shm
--stage <image>:VARUS_SIF --stage <checkout>:VARUS_SRC`. Each job copies the
genome and HISAT2 index to `/dev/shm` and the container image, code and
runlist to node-local `/tmp`, works there, and an EXIT trap copies the
results back to the shared outdir and removes both local directories, also
after a failure, `scancel` or the time-limit signal. Timings therefore
measure pyVARUS and the network, not the shared file system.

Writes to brain's Ceph `/home` and `/projects` ran at 0.4–7 MB/s, so the
copy-back leaves out everything regenerable (`.sra` files, Logan contigs,
the minimap2 index, per-run `*.tiles.npz` caches, batch and merge
directories). `VARUS.bam` is copied back only for the arms named in
`--copy-bam-arms` (here A5_logan_500); every arm keeps a `VARUS.flagstat`
computed on the node. In the smoke test, copying back the minimap2 index
and tile caches (480 MB) took 20 min.

The pre-flight check aborts a job (exit 2) unless the node has at least
`--need-tmp-gb` (80) free in `/tmp` and `--need-shm-gb` (4) free in `/dev/shm`.
A watchdog stops the run (exit 143, results still copied back) when free
space falls below `--min-tmp-gb` (30) or `--min-shm-gb` (2). With
`--prefetch-disk-gb 40`, a single job stays well within those limits.

## Datasets

* *Coelastrella tenuitheca* GCA_051903525.1 — 9 runs, no rejections (tests
  the in-loop speed-ups in isolation; Logan can prune nothing here).
* *Chlorella sorokiniana* GCA_025917655.1 (38.8 Mbp, 15 sequences) — 391 runs, 35 % of batches rejected
  in production (tests the Logan gate).
* *Drosophila melanogaster* GCF_000001215.4 — 115 457 runs, RefSeq annotation
  (tests the many-run case and intron Sn/Sp against a real annotation, as in
  the VARUS paper).
* *Mus musculus* GCF_000001635.27 (GRCm39) — 2 001 406 runs, RefSeq annotation
  (2.7 Gb genome; the runlist size forced the lazy batch order and the fresh
  pool; no A0, only the new code).

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

## Results (brain, 2026-09-23)

All 16 arms finished with exit 0. Up to 8 ran at the same time on different
nodes (snowball, batch, pinky), so they shared brain's internet uplink.
"production" is the BRAKER4 run the baseline log comes from: old code,
2 threads, final merge included. Wall time is `varus logan` plus `varus run`.
Coverage and introns are compared with A0, which uses the same code, seed
and runlist as the other arms. "Jaccard ≥ 5" uses introns with at least
5 supporting reads.

*Chlorella sorokiniana* (391 runs):

| arm | batches | rejected | wall | of which Logan | S / A0 | tiles ≥ 10 | Jaccard ≥ 5 vs A0 | vs production | vs A0 |
|---|---|---|---|---|---|---|---|---|---|
| production (2 threads) | 1000 | 343 | 4.27 h | – | 100.1 % | 7583 | 0.795 | 1.0× | 1.04× |
| A0_baseline_flags | 1000 | 371 | 4.46 h | – | 100.0 % | 7565 | 1.000 | 0.96× | 1.0× |
| A1_inloop | 1000 | 371 | 4.45 h | – | 100.0 % | 7565 | 1.000 | 0.96× | 1.0× |
| A2_parallel3 | 1000 | 387 | 1.27 h | – | 99.6 % | 7565 | 0.989 | 3.4× | 3.5× |
| A3_parallel3_prefetch | 1000 | 384 | 1.84 h | – | 99.7 % | 7565 | 0.989 | 2.3× | 2.4× |
| A3t2_parallel3_prefetch | 1000 | 381 | 1.91 h | – | 99.8 % | 7567 | 0.989 | 2.2× | 2.3× |
| A4_logan_1000 | 1000 | 246 | 1.55 h | 15 min | 102.5 % | 7589 | 0.767 | 2.8× | 2.9× |
| A5_logan_500 | 500 | 187 | 1.65 h | 17 min | 89.9 % | 7501 | 0.772 | 2.6× | 2.7× |
| A6_logan_profit | 1000 | 247 | 1.78 h | 17 min | 102.5 % | 7588 | 0.767 | 2.4× | 2.5× |

*Coelastrella tenuitheca* (9 runs):

| arm | batches | rejected | wall | of which Logan | S / A0 | tiles ≥ 10 | Jaccard ≥ 5 vs A0 | vs production | vs A0 |
|---|---|---|---|---|---|---|---|---|---|
| production (2 threads) | 1000 | 0 | 8.62 h | – | 100.0 % | 15375 | 0.766 | 1.0× | 0.53× |
| A0_baseline_flags | 1000 | 0 | 4.54 h | – | 100.0 % | 15368 | 1.000 | 1.9× | 1.0× |
| A1_inloop | 1000 | 0 | 4.49 h | – | 100.0 % | 15368 | 1.000 | 1.9× | 1.0× |
| A2_parallel3 | 1000 | 0 | 1.60 h | – | 100.0 % | 15368 | 0.994 | 5.4× | 2.8× |
| A3_parallel3_prefetch | 1000 | 0 | 0.65 h | – | 100.0 % | 15371 | 0.994 | 13.2× | 7.0× |
| A3t2_parallel3_prefetch | 1000 | 0 | 1.21 h | – | 100.0 % | 15371 | 0.993 | 7.1× | 3.8× |
| A4_logan_1000 | 1000 | 0 | 0.70 h | 1 min | 99.8 % | 15365 | 0.759 | 12.3× | 6.5× |
| A5_logan_500 | 500 | 0 | 0.40 h | 1 min | 89.0 % | 15074 | 0.724 | 21.7× | 11.4× |
| A6_logan_profit | 1000 | 0 | 0.64 h | 1 min | 99.8 % | 15365 | 0.759 | 13.4× | 7.1× |

### What the numbers show

* **Replicate noise sets the intron floor.** Production and A0 run the same
  algorithm on different random read ranges. They share only 50–59 % of all
  introns and 77–80 % of introns with ≥ 5 reads, at equal score. The
  Logan arms (0.76–0.77 at ≥ 5) sit on that floor. So the planned rule
  "≥ 90 % intron Jaccard" cannot be met, not even by a replicate.
  Arms that keep the pick order (A1–A3) reach 0.99.
* **The in-loop fixes (P1) work, but show only at low thread counts.**
  Compared with production at the same 2 threads (A3t2 vs production,
  *C. tenuitheca*), the per-batch BAM scan fell from 1.89 h to 10 min and the
  HISAT2 time from 1.74 h to 47 min. At 48 threads A1 equals A0, because
  alignment and scan are already small next to the downloads.
* **Downloads dominate everything else.** A0 spends 3.7–3.9 h of its
  4.5 h in `fastq-dump`. Three parallel downloads (A2) give 2.8–3.5×.
  Remote range dumps also get slower deep into a large run: *C. tenuitheca*
  A2's last 100 batches took 30 s median.
* **Prefetch pays off when few runs are sampled many times, and costs time
  when many runs are sampled once.**
  * *C. tenuitheca*: local `.sra` dumps take 1.3 s, so prefetch adds another
    2.5× on top of A2 (A3: 7.0× vs A0, 13.2× vs production).
  * *C. sorokiniana*: A3 is slower than A2 (1.84 h vs 1.27 h). It prefetched
    93 GB, 18 GB of it for runs whose first batch then failed the
    quality gate, and the prefetches compete with the remote dumps for
    bandwidth.
* **Logan did not speed up *C. sorokiniana*, and the benchmark shows why.**
  Only 17 of the 391 runs are usable (SRR37043103–123, 79–86 % unique
  HISAT2 alignments). All 17 are newer than the last Logan rebuild, so
  they are "absent". Nearly all other runs map below 5 % and are rejected
  after one batch each, which is 37 % of all batches.
  * Logan accepted 321 of these runs. Their contigs still cover the
    genome under `minimap2 -x splice` (2 300–4 600 tiles, yield 10–60 %).
  * Contig alignment identity shows they come from other species: the
    median `de` divergence of the contigs is **0.13–0.17**, with 0 % of
    aligned bases within 2 %. *C. tenuitheca*'s own runs score 0.0000 (100 %
    of bases within 1 %).
  * The breadth gate therefore ranked other-species runs first, and they
    took the bootstrap picks.
  * Fix: `varus logan --max-divergence` (default 0.05, see below).
* ***C. tenuitheca* has nothing for Logan to remove.** All 9 runs are good, and
  the pre-screen takes 50 s. A4 and A3 are equal within noise (0.70 h vs
  0.65 h).
* **Half the batches cost ~10 % of the score.** A5 (500 batches) reaches
  89–90 % of A0's S on both species, below the 95 % rule. For *C. sorokiniana*
  it spent 187 of its 500 batches on rejected runs.
* **`--profit-condition` never fired.** A6 ran all 1000 batches on both
  species; its differences from A4 are noise.
* **Network speed changes over the day.** *C. sorokiniana* A0 downloads took
  9.1 s median against about 5.4 s in production, so A0 comes out 4 %
  slower than production despite 48 threads. Compare arms with A0
  (same day, same infrastructure) rather than with production.
* **greif14 (no SLURM, office uplink):** *C. sorokiniana* A4 took 15.3 h
  (Logan 32 min). The run is network-bound on that machine.

### Rerun with the divergence gate and with 6 downloads (2026-09-24)

A4 and A5 on *C. sorokiniana* were rerun with `varus logan --max-divergence
0.05`, the new default. A2 was also rerun with `--parallel-downloads 6`.
All three ran on the same day as each other, but not at the same time as
the original arms.

| arm | batches | rejected | wall | of which Logan | S / A0 | tiles ≥ 10 | recovers A0 introns ≥ 10 reads | vs production | vs A0 |
|---|---|---|---|---|---|---|---|---|---|
| production (replicate) | 1000 | 343 | 4.27 h | – | 100.1 % | 7583 | 0.948 | 1.0× | 1.04× |
| A0_baseline_flags | 1000 | 371 | 4.46 h | – | 100.0 % | 7565 | 1.000 | 0.96× | 1.0× |
| A2_parallel6 | 1000 | 407 | 0.87 h | – | 99.1 % | 7564 | 0.999 | 4.9× | 5.1× |
| A4_logan_1000, divergence gate | 1000 | **7** | 1.08 h | 19 min | **105.0 %** | 7548 | 0.941 | 4.0× | 4.1× |
| A5_logan_500, divergence gate | 500 | **7** | **0.73 h** | 17 min | **95.2 %** | 7484 | 0.925 | 5.8× | **6.1×** |

* **The gate removes the wasted batches.** 321 of the 330 runs Logan
  screened were rejected as divergent, 6 had too few contigs, and 1 was
  accepted. Rejected batches in the loop fell from 246 to 7. The loop then
  samples only the 61 runs Logan could not screen, which include the 17
  usable ones.
* **More batches pass, so the score rises.** A4 with the gate reaches
  105 % of A0's score in a quarter of A0's time. A0 used only 629 of its
  1000 batches.
* **With the gate, A5 meets the old rules on *C. sorokiniana*:** 95.2 % of A0's
  score in 16 % of its wall time. It recovers 92.5 % of A0's
  well-supported introns, against 94.8 % for a replicate. On *C. tenuitheca*,
  where no batch is ever rejected, halving the batches still costs 11 %
  of the score, so the rule is met only where Logan removes waste.
* **The Logan stage is now the largest single cost of A5** (17 of 44 min).
  Almost all of that time goes into aligning contigs of runs that are then
  rejected.
* **6 parallel downloads beat 3:** 52 min against 76 min (1.46×). The
  extra in-flight picks cost 20 more rejected batches (407 against 387)
  and 0.5 % of the score. The prefetch variant (A3 with 6) was stopped by
  the watchdog on node385, whose /tmp had only 114 GB free. It was not
  rerun: `--prefetch` was retired on 2026-09-24 (see Decisions).

### *Drosophila melanogaster* (2026-09-24)

115 457 RNA-seq runs, 48 threads, all arms on the same day. The Logan arms
used the divergence gate. A4/A5 were submitted with `--prefetch` before it
was ruled out; their `varus run` time includes the prefetches. Intron
Sn/Sp use the VARUS paper's definition (Stanke et al. 2019): predicted =
distinct introns from the spliced alignments, 32 bp–350 kb (the `bam2hints`
window); reference = the 47 911 coding introns of the RefSeq annotation.
The paper reports Sn 0.935 / Sp 0.359 for *Drosophila* from 758 runs.

| arm | batches | rejected | wall | of which Logan | S / A0 | tiles ≥ 10 | runs sampled | Sn / Sp (all) | Sn / Sp (≥ 2 reads) | vs A0 |
|---|---|---|---|---|---|---|---|---|---|---|
| A0_baseline_flags | 1000 | 83 | 4.33 h | – | 100.0 % | 28979 | 265 | 0.950 / 0.199 | 0.933 / 0.293 | 1.0× |
| A2_parallel3 | 1000 | 86 | 1.58 h | – | 83.6 % | 26963 | 310 | 0.957 / 0.223 | 0.943 / 0.342 | 2.7× |
| A3_parallel3_prefetch | 1000 | – | loop 1.50 h, then hung | – | – | – | – | – | – | – |
| A4_logan_1000 (+prefetch) | 1000 | **27** | 3.94 h | 23 min | 92.8 % | 28123 | 198 | 0.957 / 0.184 | 0.936 / 0.301 | 1.0× |
| A5_logan_500 (+prefetch) | 500 | **5** | 2.65 h | 23 min | 83.6 % | 27647 | 118 | 0.930 / 0.298 | 0.896 / 0.446 | 1.4× |
| A7_defaults (6 downloads) | 1000 | 76 | **0.97 h** | – | 86.1 % | 27547 | 233 | 0.951 / 0.245 | 0.934 / **0.378** | **4.5×** |
| A8_logan_defaults (6 downloads, no prefetch) | 1000 | **18** | 2.36 h | 15.7 min | 92.7 % | 28091 | 165 | 0.954 / 0.186 | 0.932 / 0.305 | 1.8× |
| A9_defaults_sparse (A7 + sparse estimator) | 1000 | 83 | **37.6 min** | – | 74.8 % | 22622 | 260 | 0.953 / 0.265 | 0.938 / **0.410** | **6.9×** |
| A10_logan_sparse (A8 + sparse estimator, 8 connections) | 1000 | 23 | 1.03 h | 14.2 min | 92.7 % | 28105 | 170 | 0.955 / 0.185 | 0.933 / 0.302 | 4.2× |
| A11_logan_only (A10 + `--logan-only`) | 1000 | **5** | 55.8 min | 13.5 min | 93.2 % | 28143 | 144 | 0.952 / 0.246 | 0.929 / 0.366 | 4.7× |
| A12_logan_prior (A10 + uncapped tiles, prior ×40 until first batch) | 1000 | 36 | 59.0 min | 14.2 min | 84.5 % | 26109 | 161 | 0.956 / 0.190 | 0.936 / 0.323 | 4.4× |
| A9, mean of seeds 1–3 | 1000 | 71 | 39.1 min | – | 76.6 % | 23465 | 272 | 0.951 / 0.271 | 0.935 / 0.413 | 6.6× |
| A13 λ = 1, a = 0.1, mean of seeds 1–3 | 1000 | 145 | 40.1 min | – | 88.2 % | 27175 | 658 | 0.954 / 0.228 | 0.938 / 0.356 | 6.5× |
| A14 λ = 3, a = 0.1, mean of seeds 1–3 | 1000 | 88 | 39.7 min | – | **91.7 %** | 27251 | 329 | 0.949 / 0.244 | 0.933 / 0.374 | **6.5×** |

* **Three parallel downloads again give 2.7×.** A0 spends 2.9 h of its
  4.3 h in `fastq-dump`.
* **Prefetch is the reason the Logan arms are slow, and it hangs the run.**
  A3 finished its 1000 batches in 1.50 h (A2: 1.58 h) and then sat for
  3 h in 58 prefetches that had been queued near the end of the loop, until
  the job was cancelled. A4 made 388 prefetch calls for runs that are mostly
  sampled a few times each. Its loop alone would be at A2's speed.
* **The divergence gate works on a real annotation too.** Rejected batches
  fall from 83 to 27 (A4) and 5 (A5). Logan checked 547 candidates in
  23 min: 375 accepted, 94 rejected, 31 too few contigs, 47 absent. Of that
  time, 10 min was minimap2 and 12 min the BAM scan, which ran serially in
  the main process next to the download threads. The scan now runs in
  worker processes while the next chunk aligns (`--scan-workers`, default
  2). A Logan-only rerun with the same inputs (job 8227749, node202) gave
  identical counts and a stage of **13.4 min** instead of 23.1 min
  (pipeline 738 s against 1308 s; minimap2 628 s, scan 268 s of worker
  time). The scan itself also got cheaper, since it no longer shares the
  interpreter with the eight download threads.
* **Sensitivity matches the paper in every arm (0.93–0.96).** Specificity
  is lower (0.18–0.30 against 0.36) because pyVARUS keeps every intron
  seen once; the original algorithm (A0) shows the same. Requiring two
  reads brings A2 to 0.943 / 0.342, next to the paper's 0.935 / 0.359.
  The remaining differences are the read pool (115 k runs in 2026 against
  758 in 2019), HISAT2 with a growing splice DB, and the annotation version.
* **A2's score is 16 % below A0 here, against 0.4 % on *C. sorokiniana*.** A0
  drew 244 batches from one very productive run (SRR21970089); A2 spread
  its batches over 310 runs with at most 118 from any one. With 115 k
  runs the pick paths diverge early, so whether this is a cost of the
  in-flight accounting or path luck needs an A0 replicate to decide.
* **Half the batches fail the score rule again** (A5: 83.6 %).

#### Shipped defaults, with and without Logan (A7, A8)

A7 and A8 ran on 2026-09-24 with the code as shipped: `--parallel-downloads
6` by default, and for A8 `varus logan` with 16 download connections,
25-run chunks and minimap2's SAM piped into the scanner processes. No
prefetch. A7 ran on node220, A8 on node205, and both nodes were cleaned
afterwards.

* **Without Logan, six downloads give 4.5× over A0** (58 min against
  4.33 h, and 1.6× over A2's 1.58 h). The score is 86.1 % of A0 (A2:
  83.6 %). Intron sensitivity is unchanged (0.951). Specificity is the best
  of all arms: 0.378 at ≥ 2 reads, against 0.359 in the paper.
* **With Logan the run is more accurate by score but slower.** Rejected
  batches fall from 76 to 18. The score rises to 92.7 % of A0 from 165
  runs, and more tiles reach ≥ 10 reads (28 091 against 27 547). Intron
  sensitivity is the same (0.954); specificity drops to 0.186 (0.305 at
  ≥ 2 reads), as it did in A4, because A8 reports 33 % more distinct
  introns. The seeded splice DB (131 580 junctions from contigs) is a
  likely cause, but this run does not isolate it.
* **The Logan loop is limited by the estimator, not by downloads.** Mean
  time per batch in the main thread:

  | arm | download (in parallel) | align | scan | splice DB | estimate | batch wall |
  |---|---|---|---|---|---|---|
  | A7, first 100 batches | 8.2 s | 1.0 s | 0.6 s | 0.1 s | 0.2 s | 2.1 s |
  | A7, last 100 | 7.6 s | 1.2 s | 0.6 s | 0.4 s | 2.6 s | 4.8 s |
  | A8, first 100 | 8.9 s | 1.1 s | 0.7 s | 0.1 s | 2.9 s | 5.0 s |
  | A8, last 100 | 8.0 s | 1.2 s | 0.7 s | 0.6 s | 5.7 s | 8.2 s |

  With six downloads in flight, `fastq-dump` is hidden in both arms. The
  estimator costs 5.2 s per batch in A8 (1.44 h of the 2.10 h loop) and
  1.4 s in A7. With Logan, every one of the 375 accepted runs carries its
  own prior, so each is estimated and scored separately, and the cost
  grows with the observed tiles. The estimator is now the next thing to
  speed up. In A7 it also becomes the largest per-batch cost by the end
  of the run.
* **The Logan stage got slower: 15.7 min against 13.4 min** (job 8227749,
  same inputs and identical counts: 375 accepted, 94 rejected, 31 too few
  contigs, 47 absent). The expected gain did not happen:
  * 16 connections did not raise throughput. The aggregate rate stayed at
    about 8 MB/s (8 × 1.04 MB/s before, 16 × 0.49 MB/s now). The limit is
    the total bandwidth, not the rate per connection.
  * minimap2 took longer: 740 s against 628 s for the same 500 runs, even
    with 20 index loads instead of 38. Per chunk it used 10–25 of 48
    cores (CPU time / real time from minimap2's summary). The scan is
    hidden (time the main thread waited for it: 1.0 s in total). The run
    does not separate the possible causes: `-K 20M` mini-batches,
    minimap2 blocking on the pipe while the scanner reads, or the node
    (node205 against node202).

#### Sparse estimator without Logan (A9)

A9 is A7 with the sparse estimator (job 8229161, node214, copy-back rc=0).

* **The loop is 1.55× faster: 37.6 min against 58.2 min** (6.9× over A0).
  The estimate per batch fell from 1.42 s to 0.10 s on average (last 100
  batches: 2.60 → 0.13 s), and the batch wall time from 3.37 s to 2.12 s.
  With six downloads of about 8 s each in flight, the floor is about
  1.35 s per batch.
* **The score is 74.8 % of A0 against A7's 86.1 %, with the same number of
  UMRs (36.0 M against 36.1 M).** It is not the estimator. A replay of
  80 random batches through the controller agrees with a dense
  from-scratch computation of eq. 3 on every p̂, every profit (with and
  without in-flight extras) and every lazy-greedy pick.
* **The difference is which "super run" the loop happens to find.** A0
  first sampled SRR21970089 at batch 657 and then gave it 244 of the
  remaining 343 batches (9.4 M UMRs). A7 found its sibling SRR21970091 at
  batch 951 (41 batches). A9 found neither and spread its batches over
  260 runs, with at most 14 from one run. With 115 k runs and one
  bootstrap batch per unsampled run, whether and when such a run turns up
  is luck. Single-seed S differences of this size between the *Drosophila*
  arms are therefore not conclusive. The A2 bullet above already suspected
  this.
* Intron sensitivity is unchanged (0.953). Specificity is the best of all
  arms (0.410 at ≥ 2 reads), because fewer reads from any single run mean
  fewer distinct introns.

#### Sparse estimator with Logan (A10)

A10 is A8 with the sparse estimator, and Logan back at 8 download
connections without `-K` (job 8229162, node222, copy-back rc=0).

* **The loop is 2.65× faster than A8: 47.5 min against 2.10 h.** The
  estimate per batch fell from 5.2 s to 0.19 s. With the 14.2 min Logan
  stage the job takes 1.03 h (A8: 2.36 h), 4.2× over A0.
* **The Logan stage took 14.2 min** (A8: 15.7 min; L1 with the same
  settings: 13.4 min), with identical counts (375 accepted, 94 rejected,
  31 too few contigs, 47 absent). minimap2 took 698 s against L1's 628 s
  on a different node, and the S3 rate was 1.01 MB/s per connection.
* **The score is unchanged: 92.7 %,** as in A8, with rejected batches at
  23. The run went the same way: SRR36274151 got 134 batches (A8: 132).
* **The Logan loop is still 26 % slower per batch than A9** (2.75 s
  against 2.12 s). The downloads took longer in this run (10.6 s against
  8.1 s mean), which puts the floor for six in flight at about 1.8 s.
  HISAT2 (1.19 s against 1.02 s) and the splice-DB update (0.28 s against
  0.23 s) are slower because the DB starts with the 131 580 Logan
  junctions.
* **Logan against no Logan, as shipped:** 1.03 h against 37.6 min, for
  S 92.7 % against 74.8 % and 23 against 83 rejected batches. After the A9
  finding the score gap is mostly which super run was found (see above),
  not something Logan guarantees. Logan's bootstrap did put SRR36274151,
  one of its 547 candidates, into play early: its first batch was batch 16
  (A8: 9). **But its second batch came only at batch 230 (A8: 236), and
  from then on it got most batches** (the tenth at 238). For about 215
  batches the loop had evidence that this run spreads its reads twice as
  widely as a typical run, and did not act on it. That is the smoothing
  effect from the offline test: one real batch is 35 k UMRs against the
  ~326 k pseudo-counts of λ·T·p̄ + a, so after one batch the run's p̂ is
  still 90 % the pooled profile.

#### Logan-only runs and the stronger prior (A11, A12)

Both ran with the sparse estimator, next to A10 (jobs 8229241 on node267
and 8229272 on node280, copy-back rc=0). The Logan stage gave identical
counts in both (13.5 and 14.2 min).

* **A11 (`--logan-only`) is the best Logan configuration so far:** 55.8 min
  in total (loop 42.3 min), S 93.2 %, and only **5 rejected batches**
  (A10: 23, A9: 83). Dropping the runs Logan never screened costs no
  score here, because the 375 accepted runs contain SRR36274151 (152
  batches). It is still single-seed. Its second batch again came late
  (batch 222).
* A11 reports 25 % fewer distinct introns than A10 (185 k against 248 k)
  at the same sensitivity, so specificity is higher (0.366 against 0.302
  at ≥ 2 reads). This run does not show why.
* **A12 did worse: S 84.5 %.** It sampled SRR36274151 once, at batch 15,
  and never again, and spread its batches over 161 runs (the most-sampled
  had 37). The strong uncapped prior makes an unsampled Logan run's p̂
  follow its expression profile, so its predicted gain is realistic. A run
  with one real batch goes back to 90 % pooled profile under λ = 10, and
  that broad profile overstates its gain. The two kinds of run are
  therefore scored on different scales. The prior strength only makes
  sense together with the smoothing, so A12 is repeated with λ = 3,
  a = 0.1 as A16.

#### Does the Logan prior rank unsampled runs? (offline, 24 runs)

A8 gave 132 of its 1000 batches to one run (SRR36274151) and 393 to Logan's
top 50. To see what the prior actually knows, job 8229233 (node234, 4 min)
took 24 accepted *Drosophila* runs: Logan's top 8, the 8 runs A8 sampled
most, and 8 random accepted runs. For each run it aligned the Logan
contigs as `varus logan` does, plus two real 50 000-spot batches (at 25 %
and 75 % of the run, HISAT2, UMRs per 5-kb tile as in the loop).

For each run r and each contig weighting, the estimator's distribution for
an unsampled run was rebuilt, p̂_r ∝ a + λ·T·p̄ + s·q_r, with q_r the
normalised contig tile weights and s the prior mass in pseudo-UMRs. From it
came the predicted gain of one batch, Σ log1p(x + E·p̂_r) − log1p(x), with x
the A0 end coverage × 0.05, × 0.3 or × 1. This was compared with the realised
gain of the run's real batch 2, Σ log1p(x + c) − log1p(x). The number that
matters for picking is the rank correlation ρ across the 24 runs.

| predictor (x = 0.3 × A0) | ρ, E from Logan yield | ρ, true E | predicted / real gain |
|---|---|---|---|
| no Logan shape (shared prior) | 0.08 | – | – |
| **current: ka ≤ 50 × len, s = 25 k** | **0.09** | 0.02 | 3.4 |
| same weights, s = 10⁶ | 0.30 | 0.59 | 2.8 |
| uncapped ka × len, s = 10⁶ | 0.58 | 0.76 | 2.4 |
| uncapped ka, s = 10⁶ | 0.59 | 0.69 | 2.1 |
| uncapped ka, Logan shape only | 0.63 | 0.66 | 1.7 |
| the run's own batch 1, s = 35 k (loop after one download) | – | 0.09 | 3.1 |
| the run's own batch 1, shape only | – | 0.99 | 1.0 |

The pattern is the same for x = 0.05 and x = 1.

* **The shipped prior does not rank unsampled runs.** The smoothing term
  a + λ·T·p̄ adds up to T + 10·T ≈ 326 k pseudo-UMRs over *Drosophila*'s
  29.6 k tiles. The prior adds 25 k (`--logan-prior-batches 1`). So p̂ for
  an unsampled run is 93 % the pooled profile, and its predicted gain
  hardly depends on the run: ρ 0.09, no better than no shape at all.
  Capping the k-mer abundance at 50 (`--ka-cap`) also flattens the
  expression profile a read batch follows. With the cap, the predicted
  number of tiles one batch hits is 1.09 × the observed; uncapped it is
  0.87, about what the run's own batch 1 predicts (0.84).
* **Uncapped abundance plus a dominant prior brings ρ to about 0.6.** With
  the true batch size it reaches 0.7–0.8. Logan's yield predicts batch size
  well (ρ 0.71 against UMRs) but not breadth (0.09). Breadth has to come
  from the contig shape.
* **The same smoothing also swamps real data.** A downloaded run's own
  batch at its real weight (35 k UMRs) predicts its next batch with ρ 0.09.
  The same counts used as a shape give 0.99. Batches of one run are nearly
  interchangeable, but the paper's λ = 10 keeps p̂ close to the pooled
  profile until a run has about ten batches. That is the original
  algorithm, and this benchmark does not change it.
* **Coordinate-sorted runs defeat range batches.** Logan's rank-1 run
  SRR7866341 (94 % unique) put its 25 % batch into 82 tiles and its 75 %
  batch entirely into the mitochondrial genome (NC_024511.2, tile 2). The
  run's spots are stored in genome order, so a 50 000-spot range covers one
  locus, while its contigs cover 16 037 tiles. Neither the prior nor the
  first batch can tell; the loop learns it after one batch (A8 sampled it
  once).
* **SRR36274151 is a legitimate favourite.** A8 gave it 132 batches. Its
  batches have the lowest unique rate of the set (56 %), but they hit the
  most tiles (14 142 against 1 800–9 700 for the others).

**Would Logan have found A0's super runs?** None of the four runs that
the no-Logan arms exploited (SRR21970089/91, SRR29129684/85) was among the
547 screened candidates. Job 8229396 collected them the same way. Among
the 28 runs, ranked by the gain of one real batch at x = 0.3 × A0:

| run | real | predicted, uncapped ka | predicted, ka ≤ 50 × len | tiles per batch | contig tiles |
|---|---|---|---|---|---|
| SRR21970091 (A7: 41 batches) | #1 | #1 | #1 | 17 649 | 17 020 |
| SRR36274151 (A8: 132) | #2 | #6 | #16 | 14 142 | 26 317 |
| SRR21970089 (A0: 244) | #3 | #2 | #3 | 16 947 | 15 491 |
| SRR29129685 (A0: 115) | #4 | #3 | #2 | 11 181 | 25 227 |
| SRR29129684 (A7: 16) | #10 | #8 | #5 | 7 854 | 22 780 |

Pure Logan shape with E from the Logan yield; ρ over all 28 runs is 0.75
uncapped and 0.57 capped (0.73–0.81 and 0.55–0.64 over the three
backgrounds).

* **The contigs identify super runs.** The best runs are not the ones
  whose contigs cover the most tiles (15–17 k, against 22–26 k for
  others). They are the runs whose expression is even, so a batch of
  50 000 spots spreads widely. The uncapped abundance profile sees this.
* **So with a huge runlist, screening more candidates can pay off.** The
  cost is the limit: about 1.7 s per candidate, most of it the ~8 MB/s
  S3 download. A cheaper pre-rank looks possible. The effective number of
  contigs, exp(entropy of ka × length), predicts tiles per batch with
  ρ 0.72 without using genome positions. Computed on the first 10 % of a
  contig file (as an HTTP range request would return), it still gives
  ρ 0.61 (5 %: 0.56). This is from mapped contigs only; a real pre-screen
  would also need the gate for foreign runs.

**Smoothing (λ, a) decides between exploring and exploiting.** On the
28 runs, the estimator was given each run's real batch 1 as its only
observation. Predicted gain of the next batch vs the realised gain of
batch 2 (x = 0.05 × A0; x = 0.3 × A0 behaves the same):

| λ | a | share of batch 1 in p̂ | ρ | runs beating an unsampled run | rank of SRR21970091 / SRR36274151 |
|---|---|---|---|---|---|
| 10 | 1 (shipped, paper) | 10 % | 0.11 | 9 / 28 | #19 / #25 |
| 3 | 1 | 23 % | 0.44 | 0 / 28 | #3 / #21 |
| 1 | 1 | 37 % | 0.62 | 0 / 28 | #2 / #16 |
| 3 | 0.1 | 28 % | 0.75 | 2 / 28 | #2 / #7 |
| 1 | 0.1 | 52 % | 0.94 | 3 / 28 | #1 / #3 |
| 0.3 | 0.1 | 75 % | 0.99 | 0 / 28 | #1 / #3 |

* With the shipped values, the runs whose second batch the loop prefers
  over a fresh run are unrelated to the runs that are actually best
  (ρ 0.11). That explains the 215 batches A10 waited before it returned to
  SRR36274151.
* Lowering λ alone does not help. The pseudo-count a = 1 per tile
  (29.6 k in total) spreads an unsampled run's p̂ over every tile, so its
  predicted gain beats every sampled run and the loop would only explore.
* λ = 1 with a = 0.1 is the balance point: batch 1 carries half the
  weight, the next batch is predicted with ρ 0.94, and the only runs that
  beat a fresh run are the best 2–3, super runs included. A13/A14 test
  this with three seeds each against three seeds of A9.

A12 tests the change in the loop: tile weights uncapped (`varus logan
--tile-weight ka --tile-ka-cap 0`; intron weights keep the cap), a prior
of 40 batch equivalents (10⁶ pseudo-UMRs), and the prior dropped once a run
has its first batch (`--logan-prior-first-only`), so it ranks unsampled
runs without drowning the real counts afterwards.

#### Smoothing in the loop (A13, A14; three seeds each)

Jobs 8229480–87 (nodes 214, 222, 236, 265, 315). The first four (A9_s2,
A13_s1, A13_s2, A14_s1) ran before the one-pass scan was synced, the other
four after it. That changes the scan time (about 10 min against 7.3–8.7 min
per 1000 batches), not what is sampled.

| arm | seed | wall | rejected | S / A0 | tiles ≥ 10 | runs sampled | Sn / Sp (≥ 2 reads) | most-used run: batches (1st, 2nd) |
|---|---|---|---|---|---|---|---|---|
| A9 λ = 10, a = 1 | 1 | 37.6 min | 83 | 74.8 % | 22622 | 260 | 0.938 / 0.410 | SRR29129684: 14 (417, 423) |
| | 2 | 40.7 min | 58 | 73.9 % | 22043 | 303 | 0.941 / 0.396 | SRR12634507: 12 (611, 617) |
| | 3 | 39.1 min | 71 | 81.0 % | 25730 | 254 | 0.925 / 0.433 | SRR1297296: 48 (169, 193) |
| A13 λ = 1, a = 0.1 | 1 | 40.6 min | 141 | 92.1 % | 27918 | 593 | 0.941 / 0.300 | SRR35544217: 125 (780, 786) |
| | 2 | 41.9 min | 143 | 88.4 % | 27087 | 686 | 0.941 / 0.383 | SRR8645649: 89 (222, 228) |
| | 3 | 37.9 min | 152 | 84.0 % | 26521 | 696 | 0.932 / 0.385 | SRR15347356: 41 (210, 215) |
| A14 λ = 3, a = 0.1 | 1 | 40.2 min | 81 | 96.8 % | 27793 | 285 | 0.931 / 0.317 | SRR29755228: 263 (513, 518) |
| | 2 | 40.2 min | 86 | 87.7 % | 27132 | 381 | 0.944 / 0.364 | SRR29755174: 47 (754, 760) |
| | 3 | 38.6 min | 96 | 90.7 % | 26827 | 321 | 0.923 / 0.441 | SRR15347363: 114 (455, 460) |

For comparison: A0 took 244 batches from SRR21970089 (first at 657), A11
took 152 from SRR36274151 (first at 16, second at 222).

* **λ = 3, a = 0.1 raises S from 76.6 % to 91.7 % of A0 (seed means) at
  the same wall time** (39.7 against 39.1 min, 6.5× faster than A0). Every
  A14 seed beats every A9 seed, and seed 1 comes within 3 % of A0. Tiles
  with ≥ 10 UMRs go from 23.5 k to 27.3 k (A0: 29.0 k).
* **The loop now stays on good runs.** Under λ = 10 no run got more than
  48 batches; under λ = 3 the most-used run got 47–263. The gap between a
  run's first and second batch (5–6 batches) is the same in all arms: it
  is the six downloads that were already in flight when the first batch
  was scored. What changes is how long the loop keeps coming back.
* **λ = 1 explores too much.** It samples 593–696 runs, twice as many as
  λ = 3, and 141–152 batches fail the quality gate (λ = 3: 81–96, A9:
  58–83). Its S is 3.5 points lower on average.
* **Sensitivity is unchanged** (0.92–0.94 at ≥ 2 reads). Specificity at
  ≥ 2 reads falls from 0.41 to 0.37 (λ = 3) because more introns are
  predicted (122 k against 109 k at ≥ 2 reads, seed means). A0 follows the
  same pattern: the highest S, 153 k introns and Sp 0.29.
* **The seed still matters** (λ = 3: 87.7–96.8 %). Without Logan the
  run the loop ends up exploiting is found at batch 455–780, so which one
  it is, and how many batches remain for it, is luck. A15 combines λ = 3,
  a = 0.1 with Logan, which found SRR36274151 at batch 16 in A11; A16 adds
  A12's strong prior, which failed under λ = 10 because a sampled run's p̂
  went back to the pooled profile while an unsampled run's did not.

#### Bigger batches from the same run (job 8229443, node214)

Would one large batch be cheaper than several 50 000-spot batches from a
run the loop picks repeatedly? For three runs, spot ranges of 50 k, 100 k,
250 k and 500 k were each timed twice, with the sizes in opposite order
in the two rounds. The steps were `fastq-dump`, HISAT2 with 48 threads and
A9's final splice DB (129 k junctions), and the loop's scan
(`count_bam_stats` + `extract_introns_from_bam`). Seconds, mean of the two
repeats:

| run | size | download | align | scan | per 50 k: download / align / scan |
|---|---|---|---|---|---|
| SRR21970091 (paired, super run) | 50 k | 8.4 | 1.37 | 0.58 | 8.4 / 1.37 / 0.58 |
| | 500 k | 25.1 | 7.59 | 4.95 | 2.5 / 0.76 / 0.50 |
| SRR36274151 (single, super run) | 50 k | 27.1 | 0.93 | 0.44 | 27.1 / 0.93 / 0.44 |
| | 500 k | 34.1 | 4.44 | 4.62 | 3.4 / 0.44 / 0.46 |
| SRR1804047 (paired, typical) | 50 k | 5.9 | 1.13 | 0.49 | 5.9 / 1.13 / 0.49 |
| | 500 k | 14.7 | 7.07 | 5.04 | 1.5 / 0.71 / 0.50 |

(100 k and 250 k fall on the same lines.)

* **Download: 3–8× cheaper per read in large ranges.** Each `fastq-dump`
  call has a fixed cost of 5–7 s (26 s for SRR36274151, the run A8/A10
  took 130+ batches from), plus 1.5–4 s per 100 k spots.
* **HISAT2: about half of a 50 k batch is fixed cost** (~0.5 s at 48
  threads), so 500 k spots cost 0.44–0.76 s per 50 k instead of
  0.93–1.37 s.
* **The scan is linear** (0.45–0.5 s per 50 k at every size) and gains
  nothing. It read each BAM twice (UMR counts, then introns). Since
  2026-09-24 it is one pass (`tiles.scan_batch_bam`): on a 100 k-record
  slice of a *Drosophila* `VARUS.bam`, 0.28 s against 0.40 s with a warm
  cache, with identical UMR, spliced-read and intron counts.
* **Effect on the loop.** With six downloads in flight, the download is
  already hidden, so the gain is HISAT2's fixed cost plus the per-batch
  splice-DB update and estimate (~0.35 s), about 0.8 s per merged 50 k.
  If the ~60 % of A9's batches that went to repeatedly picked runs came in
  250 k batches, the loop would drop from 37.6 to about 30 min (−20 %).
  Runs with a slow first byte like SRR36274151 gain the most, and with
  Logan they are exactly the runs the loop exploits.
* **Overlapping HISAT2 for the next batch with the scan of the current
  one** saves about as much for every batch (main thread ~2.0 → ~1.2 s)
  without changing what is sampled. Its limit is then the download
  throughput (six streams, ~1.35–1.8 s per batch), which merged batches
  relieve. The two are complementary. Implemented on 2026-09-25 as
  align-ahead: −23 % (A17/A18 below).

#### Logan with λ = 3 (A15, A16; 2026-09-25)

A15 is A11 (`--logan-only`) with λ = 3, a = 0.1, three seeds. A16 adds
A12's uncapped tile weights and a ×40 prior until a run's first batch.
All four ran before align-ahead.

| arm | Logan stage | loop | s/batch | download per call | rejected | batches from SRR36274151 (first at) | S / A0 | Sn / Sp | Sn / Sp (≥ 5) |
|---|---|---|---|---|---|---|---|---|---|
| A11 (λ = 10) | 13.5 min | 41.2 min | 2.47 | 10.1 s | 5 | 152 (16) | 93.2 % | 0.952 / 0.246 | 0.873 / 0.531 |
| A15_s1 | 13.3 min | 57.6 min | 3.46 | 17.6 s | 6 | 518 (17) | 102.2 % | 0.943 / 0.235 | 0.852 / 0.442 |
| A15_s2 | 14.3 min | 56.8 min | 3.41 | 17.4 s | 4 | 512 (16) | 102.1 % | 0.943 / 0.237 | 0.854 / 0.444 |
| A15_s3 | 13.6 min | 57.3 min | 3.44 | 17.5 s | 3 | 510 (16) | 102.1 % | 0.943 / 0.236 | 0.853 / 0.444 |
| A16 | 13.4 min | 56.2 min | 3.37 | 17.5 s | 2 | 527 (19) | 104.1 % | 0.928 / 0.239 | 0.812 / 0.429 |

* **Best score so far, and no seed spread.** A15 reaches 102 % of A0 in
  every seed (A14 without Logan: 87.7–96.8 %). Logan's bootstrap finds the
  super run SRR36274151 at batch 16–19 and λ = 3 keeps the loop on it:
  about half of all batches. Rejected batches fall to 2–6.
* **But 40 % slower than A11 and 1.8× slower than A14.** SRR36274151 has
  a 27 s fixed cost per `fastq-dump` call (see "Bigger batches"); with half
  the batches from it, the mean download per call rises from 10 to 17.5 s
  and six streams deliver a batch only every ~3.4 s. The loop is download
  bound, so align-ahead would not help much here. Total with the Logan
  stage: ~72 min, 3.6× over A0.
* **Intron sensitivity at ≥ 5 reads drops** (0.853 against 0.873–0.910):
  depth concentrates on one run's genes. S counts tiles, not introns, so
  the score rewards this; the downstream effect (GeneMark-ETP, BUSCO) is
  not measured.
* **A16's strong prior adds a second exploited run** (SRR21684949, 115
  batches), 2 points of S, and loses more intron sensitivity (0.812 at ≥ 5).
* **This is exactly the case merged batches address:** one 500 k-spot call
  to SRR36274151 took 34 s against 27 s for 50 k, i.e. 3.4 s per 50 k
  instead of 27 s. With the repeated picks of this run merged, the A15 loop
  would be limited by alignment again (~30 min), for ~45 min with Logan.

#### Align-ahead (A17 vs A18, 2026-09-25)

A17 aligns the next batch in a background thread while the main thread
scans and scores the current one; A18 is the same code with
`--no-align-ahead`. Both use A14's λ = 3, a = 0.1. All six jobs started
within a minute of each other on separate nodes. Loop times are from `BatchTimings.tsv` (seconds per
batch are means over 1000 batches).

| arm | loop | run wall | s/batch | wait | align | scan | db + est | rejected | S / A0 | speedup vs A0 |
|---|---|---|---|---|---|---|---|---|---|---|
| A17_ahead_s1 | 28.7 min | 30.9 min | 1.72 | 0.38 | 1.31 | 0.64 | 0.60 | 59 | 99.3 % | 8.4× |
| A17_ahead_s2 | 29.5 min | 32.0 min | 1.77 | 0.48 | 1.32 | 0.72 | 0.48 | 57 | 98.3 % | 8.1× |
| A17_ahead_s3 | 27.6 min | 30.1 min | 1.66 | 0.59 | 1.11 | 0.60 | 0.39 | 90 | 90.9 % | 8.6× |
| A18_noahead_s1 | 36.6 min | 38.9 min | 2.20 | 1.21 | 1.17 | 0.50 | 0.40 | 86 | 96.5 % | 6.7× |
| A18_noahead_s2 | 39.6 min | 42.0 min | 2.38 | 1.33 | 1.27 | 0.54 | 0.42 | 59 | 96.8 % | 6.2× |
| A18_noahead_s3 | 35.1 min | 37.5 min | 2.10 | 1.17 | 1.09 | 0.51 | 0.35 | 88 | 91.2 % | 6.9× |

* **−23 % loop time for the same seed on the same code** (s1 36.6 → 28.7,
  s2 39.6 → 29.5, s3 35.1 → 27.6 min). Three-seed mean run wall 31.0 min
  (A17) against 39.5 min (A18): 8.4× over A0 instead of 6.6×.
* **Picks are not affected.** Seed 3 took essentially the same path in both
  arms (S 90.9 vs 91.2 %, 149 k introns in both). Seeds 1 and 2 differ in S
  because pick timing shifts and which super run is found is luck
  (A14_s2 87.7 %, A18_s2 96.8 %, A17_s2 98.3 %), not because of align-ahead. A18
  reproduces A14 (s1: 38.9 vs 40.2 min, S 96.5 vs 96.8 %).
* **Less than the 2.0 → 1.2 s expected** (measured 2.23 → 1.72 s). HISAT2
  (1.1–1.3 s, now in the background) is as long as the main thread's
  work (scan + DB + estimate, 1.1–1.3 s). Both slow down by ~20 % because
  they now share the 48 cores (HISAT2 `-p 48` plus `samtools sort -@ 47`).
  The remaining 0.4–0.6 s of wait is the main thread waiting for the
  aligner or for a download.
* Sn/Sp against RefSeq move with the picks, not with the arm (s3:
  0.943/0.308 vs 0.941/0.308).

#### Thread budget (A19 vs A17, 2026-09-25)

A19 is A17 with the thread budget (HISAT2 `-p 46` instead of 48,
`samtools sort -@ 4` instead of 47, rolling merge `-@ 4` instead of 24 and
taken out of the aligner's share while it runs). Same seeds, run 50 min
after A17.

| seed | loop A17 → A19 | align (s/batch) | scan (s/batch) | wait (s/batch) |
|---|---|---|---|---|
| 1 | 28.7 → 28.4 min | 1.31 → 1.27 | 0.64 → 0.62 | 0.38 → 0.41 |
| 2 | 29.5 → 29.8 min | 1.32 → 1.23 | 0.72 → 0.63 | 0.48 → 0.61 |
| 3 | 27.6 → 27.4 min | 1.11 → 1.06 | 0.60 → 0.58 | 0.59 → 0.61 |

* **Time-neutral**: the loop is within ±1 % per seed. Alignment and scan
  are 3–9 % faster each, but the gain goes into waiting for downloads.
* So the ~20 % slowdown of HISAT2 and the scan under align-ahead was not
  thread oversubscription. A 48-CPU allocation on brain is 24 physical cores
  with two hardware threads each (`ThreadsPerCore=2`), so the scan shares a
  core with HISAT2 threads either way.
* Kept: it costs nothing and stops concurrent stages from starving each
  other at lower thread counts or with Logan's scanner processes.

#### Merged batches (A20, A21 vs A22, A19; 2026-09-25)

`--merge-batches 10`: after a run's first batch passed the quality gate,
a repeated pick claims the following spot ranges (up to 10 × 50 k spots,
one `fastq-dump` call, one alignment, one scan) while greedy selection
would pick the same run again with the claimed batches counted as
observed. All arms on the current code (align-ahead, thread budget),
λ = 3, a = 0.1, submitted together.

| arm | Logan stage | loop | total (sacct) | downloads (merged) | runs sampled | rejected | S / A0 | Sn / Sp ≥ 5 | speedup vs A0 |
|---|---|---|---|---|---|---|---|---|---|
| A22 Logan, no merge | 13.3 min | 54.5 min | 69.6 min | 1000 (0) | 157 | 4 | 102.3 % | 0.851 / 0.441 | 3.8× |
| **A20_s1 Logan + merge** | 13.3 min | 19.5 min | 34.2 min | 362 (114) | 127 | 7 | 100.7 % | 0.869 / 0.452 | 7.6× |
| A20_s2 | 14.4 min | 18.5 min | 34.3 min | 336 (114) | 124 | 3 | 100.9 % | 0.875 / 0.450 | 7.6× |
| A20_s3 | 13.5 min | 19.6 min | 34.5 min | 365 (122) | 127 | 6 | 101.1 % | 0.867 / 0.453 | 7.6× |
| A19_s1 no Logan, no merge | – | 28.4 min | 30.9 min | 1000 (0) | 243 | 56 | 98.0 % | 0.898 / 0.400 | 8.5× |
| A19_s2 | – | 29.8 min | 32.8 min | 1000 (0) | 448 | 77 | 83.7 % | 0.912 / 0.531 | 8.0× |
| A19_s3 | – | 27.4 min | 30.1 min | 1000 (0) | 383 | 81 | 91.7 % | 0.872 / 0.586 | 8.7× |
| **A21_s1 merge** | – | 18.9 min | 20.6 min | 257 (101) | 129 | 27 | 85.4 % | 0.903 / 0.541 | 12.8× |
| A21_s2 | – | 17.3 min | 19.1 min | 230 (105) | 111 | 23 | 99.0 % | 0.891 / 0.515 | 13.8× |
| A21_s3 | – | 16.7 min | 18.6 min | 212 (104) | 97 | 19 | 90.5 % | 0.876 / 0.582 | 14.4× |

* **With Logan the loop is 2.8× faster** (54.5 → 18.5–19.6 min) and the
  whole run takes 34 min instead of 70 (7.6× over A0). A22's main thread
  waited 38 of its 55 min for downloads of SRR36274151 (27 s per call);
  merged, 10 batches of it cost ~25 s. S drops 1.4 points (102.3 →
  100.9 %), intron sensitivity at ≥ 5 reads rises (0.851 → 0.867–0.875).
* **Without Logan the loop is 37 % faster** (28.5 → 17.6 min, 13–14×
  over A0). Mean S is unchanged (91.1 → 91.6 %); the seed spread stays,
  because finding a super run is still luck without Logan.
* **Fewer runs are explored.** A21 sampled 97–129 runs against 243–448:
  1000 batches go in ~230 downloads, and each merged pick is decided
  against the estimate at pick time, before its data arrive. Rejected
  batches fall accordingly (19–27 against 56–81). Tiles ≥ 1 barely move
  (28.7 k vs 28.9 k), and intron Sn is unchanged.
* **Where the time goes now:** the main thread's BAM scan (7.5–12.6 min
  per 1000 batches, ~0.9 s per 50 k at 500 k-spot batches, single-threaded
  Python) is the largest part of the loop, followed by waiting (3–8 min).
  With Logan, its 13–14 min stage is now 40 % of the run.
* `--merge-batches 10` is the default since 2026-09-25 (not used with
  `--parallel-downloads 1`).

#### Parallel scan of merged batches (A23, A24 vs A20, A21; 2026-09-25)

`--scan-workers 4`: a merged batch's BAM is indexed (`samtools index -@ 2`)
and scanned by genome region in 4 spawned processes (16 regions on tile
boundaries); a read name whose records all lie in one region with `NH:i:1`
is finished there, all others are resolved in the main process. Single
batches keep the one-pass scan. The aligner gives up 4 cores instead of 1
for the scanners (48 threads: HISAT2 `-p 43`).

| arm | merged scan | loop | t_wait (sum) | S / A0 | Sn / Sp ≥ 5 | rejected |
|---|---|---|---|---|---|---|
| A20 Logan + merge (s1–s3) | 3.1–3.6 s | 19.6–20.7 min | 7.1–8.4 min | 100.7–101.1 % | 0.867–0.875 / 0.450–0.453 | 3–7 |
| A23 + parallel scan | 1.8–1.9 s | 19.2–20.4 min | 9.4–10.1 min | 100.7–101.2 % | 0.866–0.871 / 0.444–0.452 | 4–5 |
| A21 merge, no Logan | 5.9–6.8 s | 18.1–20.4 min | 3.2–4.2 min | 85.4–99.0 % | 0.876–0.903 / 0.515–0.582 | 19–27 |
| A24 + parallel scan | 2.0–3.5 s | 16.5–19.4 min | 6.6–9.8 min | 77.2–89.3 % | 0.879–0.914 / 0.483–0.581 | 26–45 |

* **Exact on real data.** Job 8232714 aligned three merged-size *Drosophila*
  batches (SRR18131184, SRR38033199, SRR21125149; 244–488 k paired reads)
  and scanned each with the one-pass scan and with 4, 16 and 64 regions:
  UMR tile counts, read and spliced-read counts and intron multiplicities
  were identical in all nine comparisons. Scan time 3.7 → 1.0 s (16 regions).
* **The scan is 1.8–3× faster, the loop is not.** The time saved moves
  into waiting for downloads (`t_wait` +2 to +5 min): with merged batches
  the loop is download bound on brain. Expected gain on a faster network or
  with fewer threads per download: up to the scan's 5–12 min per 1000
  batches.
* **A24's lower S is not the scan.** A24 received as many UMRs as A21
  (38.8–39.5 M vs 38.9–40.3 M) but seeds 1 and 2 spread them over 205–218
  runs instead of 97–129 and landed on fewer high-coverage tiles. The
  non-Logan arms scatter like this in every configuration (A14 88–97 %,
  A19 84–98 %, A21 85–99 %) because download completion order, and with it
  the pick sequence, depends on timing. With Logan the three seeds stay
  within 0.5 points in every arm.
* `--scan-workers 4` is the default since 2026-09-25 (0 or 1 = one-pass
  scan; only used with merged batches).

#### Logan stage with 3 minimap2 groups (A25 vs A23; 2026-09-25)

A single scanner process reading minimap2's SAM stream throttles minimap2:
one 25-run chunk took 28.4 s with minimap2 writing to `/dev/null` and 38.7 s
piped into one scanner (job 8232497). Split into 2 or 3 minimap2 processes
with a half or a third of the threads, each piped into its own scanner, the
chunk took 32.4 and 30.8 s. `varus logan --align-groups 3` splits each chunk
into 3 groups balanced by contig bp; a run is never split, and minimap2
aligns each contig independently, so the per-run statistics cannot change.

| arm | Logan stage | loop | total | `LoganRanking.tsv` | S / A0 | Sn / Sp ≥ 5 |
|---|---|---|---|---|---|---|
| A23 (1 group, s1–s3) | 14.1–14.4 min | 19.2–20.4 min | 33.7–35.2 min | – | 100.7–101.2 % | 0.866–0.871 / 0.444–0.452 |
| A25 (3 groups) | 11.4–11.6 min | 18.2–19.8 min | 30.1–31.9 min | byte-identical to A23, all 3 seeds | 100.7–101.1 % | 0.865–0.869 / 0.448–0.450 |

* **Logan stage −20 %** (855 → 688 s): minimap2 wall 697 → 501 s. Thread
  budget at 48 threads: minimap2 3 × 14, 3 scanners, 1 for main and
  downloads. `timings_s.scan` in `logan_summary.json` now sums the three
  scanners' times (1431 s), so it exceeds the wall time.
* **Rankings unchanged:** `LoganRanking.tsv` is byte-identical to A23 in all
  three seeds, as the unit test (`groups 1 vs 3`) predicts.
* **Whole run 8.3–8.7× over A0** (A0 4.33 h; A25 30–32 min including the
  Logan stage) with S 101 % of A0. The stage is now 37 % of the run;
  its remaining floor is the download of ~4.5 GB of contigs at ~8 MB/s
  (8 connections), about 9 min on brain.
* `--align-groups 3` is the default at 48 threads since 2026-09-25. Each
  minimap2 process loads its own index, so the count is capped by memory
  (see the mouse test below), not by genome size.

### Final settings on the algae: λ and Logan (B1–B4; 2026-09-25)

Jobs 8233237–8233260, 3 seeds per arm, all 24 at the same time. No
`--logan-only`: Logan could screen only 330 of *C. sorokiniana*'s 391 runs and
the rest include usable ones. S is relative to the A0 of 2026-09-23; loop
and Logan times are ranges over seeds.

| | Logan stage | loop | rejected | S / A0 | tiles ≥ 10 | introns |
|---|---|---|---|---|---|---|
| ***C. tenuitheca*** (production 8.25 h, 1000 batches) | | | | | | |
| B1 Logan, λ = 3 | 0.7 min | 14.0–14.2 min | 0 | 100.0–100.2 % | 15 384–15 395 | 193.5–193.9 k |
| B2 λ = 3 | – | 15.8–16.6 min | 0 | 100.4 % | 15 396–15 400 | 193.6–193.7 k |
| B3 Logan, λ = 10 | 0.7 min | 13.6–13.9 min | 0 | 99.7–99.8 % | 15 364–15 388 | 193.8–194.4 k |
| B4 λ = 10 | – | 13.9–14.0 min | 0 | 99.7–99.9 % | 15 351–15 364 | 194.0–194.3 k |
| ***C. sorokiniana*** (production 4.14 h, 343 rejected) | | | | | | |
| B1 Logan, λ = 3 | 8.4–8.8 min | 14.3–23.7 min | 11–13 | **105.9–106.8 %** | 7 585–7 610 | 173.8–179.8 k |
| B2 λ = 3 | – | 23.1–25.8 min | 402–420 | 99.5–99.8 % | 7 565–7 574 | 144.3–146.1 k |
| B3 Logan, λ = 10 | 8.4–8.6 min | 14.8–22.6 min | 11–13 | 105.8–106.7 % | 7 585–7 612 | 174.6–180.2 k |
| B4 λ = 10 | – | 19.2–20.2 min | 403–414 | 98.6–98.8 % | 7 554–7 564 | 143.7–144.4 k |

* **λ = 3 is never worse.** *C. tenuitheca* +0.3–0.6 points of S with and
  without Logan, *C. sorokiniana* +1 point without Logan and equal with it
  (Logan leaves 68 runs), *Drosophila* +15 points (A14). λ = 3, a = 0.1 is the
  default since 2026-09-25.
* **Logan pays on *C. sorokiniana*, not on *C. tenuitheca*.** It cuts rejected
  batches from 402–420 to 11–13 and raises S from 99 to 106 % and the intron count by
  20 %. Its stage costs 8.6 min, so the total (23–32 min) is 3–12 min
  longer than without Logan (19–26 min). On *C. tenuitheca* (9 runs, none
  foreign) it costs 0.7 min and changes S by less than the seed scatter.
* **With Logan the loop concentrates on few runs.** In B1 seed 1, 749
  downloads went to 2 runs, in seed 3 to 3 runs; without Logan (B4 seed 1)
  20 runs share them. Logan keeps 1 accepted and 67 unscreened runs, whose
  expected reads are weighted by the gate's acceptance rate. S and the
  intron count are higher, but read diversity is lower.
* **The loop is download bound.** The Logan seeds 1 and 2 took 23 min against
  14 min for seed 3 because they drew most batches from a run whose
  downloads were slow (1.6–1.8 h of summed download time against 38–40 min).
* Against the production runs: *C. tenuitheca* 8.25 h → 14–17 min (30–35×),
  *C. sorokiniana* 4.14 h → 23–32 min with Logan (8–11×), 19–26 min without it.

### Thread scaling (T4–T16 vs B1_s1, B2_s1; *C. sorokiniana*, 2026-09-25)

B1_s1 (Logan) and B2_s1 (no Logan) repeated with `--threads` 4, 8 and 16
and the automatic `--scan-workers` and `--align-groups` (0/0/2 scan
workers, one minimap2 process). Seed 1, snowball, submitted together.

| `--threads` | Logan stage | `varus run` with Logan | `varus run` without Logan | rejected (Logan / none) | S / A0 (Logan / none) | introns (Logan / none) |
|---|---|---|---|---|---|---|
| 4 | 64.7 min | 37.1 min | 32.6 min | 13 / 413 | 105.9 / 99.5 % | 175 k / 145 k |
| 8 | 35.6 min | 26.0 min | 25.7 min | 13 / 409 | 105.9 / 99.6 % | 175 k / 145 k |
| 16 | 19.4 min | 23.2 min | 24.6 min | 13 / 407 | 105.9 / 99.7 % | 175 k / 145 k |
| 48 | 8.7 min | 22.8 min | 25.4 min | 13 / 420 | 106.0 / 99.5 % | 175 k / 144 k |

* **Results do not depend on the thread count.** S, rejected batches and
  intron counts stay within the seed-to-seed range at every setting.
* **The read loop is download bound down to 8 threads** (`varus run`
  23–26 min). At
  4 threads HISAT2 shows (29 min summed against 12–16 min), so the loop
  takes 33–37 min.
* **The Logan stage is bound by minimap2's CPU.** With one minimap2
  process its alignment took 3821, 2074 and 1112 s with 3, 6 and 13
  minimap2 threads, i.e. 11.5–14.5 k CPU seconds each time. The contig
  downloads (1.6 k thread seconds over 8 connections) are not the limit.
  With few cores, Logan therefore costs more time than it does on the
  48-thread nodes. At 8 threads it adds 36 min to a 26 min run.

### *Mus musculus* (B2, B4–B6; 2026-09-25)

GRCm39 (2.7 Gb), 1 997 168 runs after the colorspace filter, 1000 batches,
48 threads, all jobs started at the same time. The arms are B2 (no Logan,
λ = 3; three seeds), B4 (no Logan, λ = 10), B5 (`--logan-only`, λ = 3; three
seeds) and B6 (`--logan-only`, λ = 10). No A0 was run, so S is relative to
B2_s1. Sn/Sp are measured against the 206 131 coding introns of the RefSeq
annotation.

| arm | Logan stage | `varus run` | job (sacct) | rejected | S / B2_s1 | introns | Sn / Sp | Sn / Sp ≥ 5 |
|---|---|---|---|---|---|---|---|---|
| B2_s1 | – | 15.1 min | 16.5 min | 6 | 100.0 % | 367 k | 0.831 / 0.467 | 0.744 / 0.762 |
| B2_s2 | – | 19.9 min | 21.3 min | 12 | 101.6 % | 418 k | 0.834 / 0.411 | 0.731 / 0.729 |
| B2_s3 | – | 17.0 min | 18.4 min | 10 | 89.9 % | 388 k | 0.858 / 0.455 | 0.771 / 0.752 |
| B4 (λ = 10) | – | 15.7 min | 17.9 min | 10 | 85.9 % | 361 k | 0.832 / 0.475 | 0.746 / 0.755 |
| B5_s1 | 75.6 min | 18.8 min | 95.4 min | 6 | 96.8 % | 484 k | 0.908 / 0.387 | 0.835 / 0.694 |
| B5_s2 | 77.4 min | 18.6 min | 97.7 min | 6 | 92.1 % | 490 k | 0.909 / 0.382 | 0.836 / 0.690 |
| B5_s3 | 73.8 min | 18.4 min | 93.9 min | 7 | 90.4 % | 478 k | 0.906 / 0.391 | 0.829 / 0.693 |
| B6 (λ = 10) | 75.6 min | 18.2 min | 95.8 min | 6 | 91.6 % | 493 k | 0.909 / 0.380 | 0.836 / 0.683 |

* **Speed.** `varus run` takes 15–20 min for 1000 batches, the same as on
  the algae and faster than on *Drosophila*. The runlist has 2 M runs, but
  thanks to the lazy batch order and the fresh pool the estimator and picks
  cost 3.0–3.6 min in total without Logan. With Logan they cost 8.5–9.1 min,
  because 259 runs carry priors over ~350 k tiles.
* **Logan raises intron sensitivity, not S.** Sn rises from 0.831–0.858
  to 0.906–0.909, and at ≥ 5 reads from 0.73–0.77 to 0.83–0.84. The loop
  also finds 480–490 k introns against 360–420 k. S is 90–97 % of B2_s1
  against 90–102 % without Logan: Logan spreads the reads over more genes
  instead of piling them on the best-covered tiles. Sp is lower (0.38–0.39
  against 0.41–0.48). Part of that is expected, because Sp counts every
  junction outside the coding introns as false, including UTR and
  non-coding junctions from the additional tissues.
* **Few rejected batches either way** (6–12). Unlike on *C. sorokiniana*, the
  mouse runlist has few foreign runs, so the Logan gate has little to
  remove. Of 500 screened runs, 259 were accepted, 228 rejected, 85 absent
  and 13 had too few contigs.
* **λ.** Without Logan, λ = 3 gives 90–102 % against 86 % for λ = 10 (one
  seed); with Logan the two are equal (90–97 % against 92 %), with the same
  intron Sn. This is consistent with the algae and *Drosophila*.
* **The Logan stage makes the job 5× longer** (94–98 against 17–21 min).
  The index build takes 66 s and the downloads take 29 min spread over 8
  connections. The stage was limited by its scanner: these jobs still used
  the old rule of one minimap2 process per chunk for genomes over 1 Gb, and
  the single scanner took as long as minimap2 (scan 4039 s, align 3978 s,
  pipeline 4096 s). See the next section for the memory-based group count.

### Mouse Logan stage with memory-based `--align-groups` (L1; 2026-09-25)

Job 8234577 (node234, 48 threads, `--mem 120G`), Logan stage only, seed 1,
same runlist and settings as B5_loganonly_lam3_s1. The >1 Gb rule was
replaced by a memory cap: each minimap2 process needs the index size plus
2 GiB within 60 % of the available memory.

| | B5_s1 (1 × 42 threads, 2 scanners) | L1 (3 × 14 threads, 3 scanners) |
|---|---|---|
| Logan stage | 4524 s (75.4 min) | 4006 s (66.8 min) |
| pipeline / align | 4096 / 3978 s | 3568 / 3444 s |
| scan (summed over scanners) | 4039 s | 9708 s |
| download thread-seconds | 14 114 | 13 864 |
| peak memory (MaxRSS) | – | 65.8 GB (index 10.8 GB) |

* **Same results.** `LoganRanking.tsv` and `Runlist.logan.tsv` are
  byte-identical to B5_s1.
* **11 % faster.** The scanners no longer limit the stage; minimap2 now
  does (pipeline ≈ align). The remaining 57 min are minimap2 CPU time at
  42 threads, which only more cores shorten.
* One seed, one node; brain's file system varies between jobs, so the
  gain is an estimate. The job was marked FAILED only because the job
  script's final `samtools flagstat` line returned 1 without a
  `VARUS.bam`; `varus logan` exited 0.

### Measured vs expected

| | expected (plan) | measured |
|---|---|---|
| P1 alone, 2 threads | 8.6 → 6.6 h | scan 1.89 h → 10 min, align 1.74 h → 47 min (*C. tenuitheca*) |
| P1 + P2 + P3, 2 threads | 3–4× (8.6 → ~2.2 h) | *C. tenuitheca* 7.1× (1.21 h); *C. sorokiniana* 2.2× (1.91 h) |
| + real thread count | ~6× | *C. tenuitheca* 13.2× (A3); *C. sorokiniana* 3.4× (A2) |
| + Logan | −35 % batches where foreign runs exist | rejections 371 → 246; good runs absent from Logan, other species accepted |
| + half the batches | ÷2 if the score holds | ÷1.8, but S 89–90 % (fails the 95 % rule) |
| 6 downloads as default (*Drosophila*, 48 threads) | faster than 3 | 4.5× over A0, 1.6× over 3 downloads; S 86 %, Sn/Sp unchanged or better |
| Logan: SAM piped into the scanner, 16 connections, 25-run chunks | stage −1.5 to −3 min | stage +2.3 min (13.4 → 15.7 min); throughput unchanged at ~8 MB/s |
| Logan + 6 downloads (*Drosophila*) | −rejected batches, faster loop | rejected 76 → 18, S 86 → 93 %, but 2.4× slower than A7 (estimator 5.2 s/batch) |
| Sparse estimator (*Drosophila*) | same picks, faster loop | identical picks in the replay; A7 58 min → A9 38 min (estimate 0.10 s/batch) |
| λ = 1, a = 0.1 (offline ρ 0.94) | loop exploits super runs | S 76.6 → 88.2 % (3 seeds), but twice the runs sampled and 145 rejected batches |
| λ = 3, a = 0.1 (offline ρ 0.75) | same, less exploration | S 76.6 → 91.7 % (3 seeds) at unchanged wall time; rejections unchanged |
| One-pass batch scan | −30 % scan time | 0.40 → 0.28 s per 100 k records; ~2 min per 1000 *Drosophila* batches |
| Logan + λ = 3 (A15, 3 seeds) | Logan's early finds kept by the smoothing | S 102 % of A0 in every seed, 3–6 rejected; loop 57 min (download bound: half the batches from a run with 27 s per call), 72 min with Logan |
| Merged batches (`--merge-batches 10`) | Logan + λ = 3 loop 57 → ~30 min | 54.5 → 19 min (34 min with Logan, 7.6× over A0, S 101 %); without Logan 28.5 → 17.6 min (13–14×), S unchanged |
| Align-ahead (next batch's HISAT2 during the scan) | main thread 2.0 → 1.2 s per batch | 2.23 → 1.72 s per batch, loop −23 % (paired, same code, 3 seeds); *Drosophila* 39.7 → 31.0 min, 8.4× over A0 |
| Parallel scan of merged batches (`--scan-workers 4`) | loop −5 to −10 min per 1000 batches | scan 3.1–6.8 → 1.8–3.5 s per merged batch, identical results; loop unchanged (19–20 min) because the saved time goes into waiting for downloads |
| Logan: 3 minimap2 groups per chunk (`--align-groups 3`) | chunk 38.7 → 30.8 s (−20 %) | Logan stage 14.3 → 11.5 min (−20 %), rankings byte-identical; whole run 30–32 min, 8.3–8.7× over A0, S 101 % |

### Decisions

* **Do not lower `--max-batches` in BRAKER4.** A5 fails the score rule on
  both species.
* **BRAKER4 wrapper defaults:**
  * `--parallel-downloads` defaults to 6 since 2026-09-24 (*C. sorokiniana*:
    1.46× over 3 for 20 more rejected batches and 0.5 % of the score), so
    the wrapper passes nothing; K=1 still reproduces the v1 pick sequence;
  * pass the real thread count to the wrapper;
  * do **not** set `--prefetch` (decided 2026-09-24). It wins only for
    species with a handful of runs (*C. tenuitheca*). With many runs it fills
    node-local disk with `.sra` files of runs that are then rejected
    (*C. sorokiniana* 93 GB; watchdog stop on node385), competes with the range
    dumps for bandwidth, and prefetches queued near the end keep
    downloading after the last batch: *Drosophila* A3 finished its 1000
    batches in 1.5 h (same as A2) and then hung for 3 h in 58 queued
    prefetches until it was cancelled. The flag was removed on 2026-09-26.
* **Use Logan with the divergence gate at 1000 batches.** On *C. sorokiniana*
  it removed 97 % of the wasted batches (4.1× faster than A0, score
  105 %). On *C. tenuitheca* it costs 50 s. Its value is removing runs from the
  wrong species, which the breadth gate alone could not do. On *Drosophila*
  it cut rejected batches from 83 to 27 and the intron sensitivity against
  the RefSeq annotation matches the VARUS paper. The stage costs 13–16 min
  there (23 min before `--scan-workers`, 11.5 min with `--align-groups 3`). With the shipped defaults,
  though, the *Drosophila* run is 2.4× slower with Logan than without (A8
  2.36 h against A7 0.97 h) because the estimator limits the loop. The
  gain is score (93 % against 86 % of A0), not time. On species without
  foreign runs, Logan is a speed-up only once the estimator is faster.
* **Keep 1000 batches.** 500 batches meets the score rule only when Logan
  removes waste (*C. sorokiniana* 95.2 %), not in general (*C. tenuitheca* 89 %).
* **Smoothing: λ = 3, a = 0.1 is the default since 2026-09-25** (A14: +15
  points of S on *Drosophila* at unchanged wall time; B1–B4: same or higher
  S on *C. tenuitheca* and *C. sorokiniana*, with and without Logan, three seeds
  each).
* **Replace the intron rule.** "≥ 90 % intron Jaccard" becomes "intron
  Jaccard at ≥ 5 reads no lower than the replicate floor (production vs
  A0)".

The raw tables come from `benchmark_varus.py report --outdir runs/<sp>
--baseline-log data/<sp>/baseline_varus.log`. The report reads
`baseline_Coverage.csv` and `baseline_introns.gff` next to the log.

## Logan contigs and unitigs as transcript evidence (*Drosophila*, 2026-09-25)

Question: can Logan assemblies replace sampled reads as transcript
evidence (spliced alignment → StringTie → ORFs), and do they carry
alternative splicing? Reference: the 58 674 introns of RefSeq mRNAs, split
into introns of single-isoform genes, constitutive introns (in every
isoform), and alternative ones; "competing" alternative introns overlap a
different intron of the same gene (skipping, alternative donor/acceptor).
Scripts and intermediate files are not in the repo (session scratchpad and
`~/varus_bench/logan_branches/dmel` on brain).

### Contigs of 375 accepted runs vs A8 reads

Introns from `logan/logan_introns.gff` (contigs) and `introns.gff` (reads)
of A8, compared by coordinates:

| class | n | contigs | reads |
|---|---|---|---|
| single-isoform gene | 18 062 | 98.0 % | 94.8 % |
| constitutive | 16 798 | 99.0 % | 99.0 % |
| alternative, non-competing | 4 054 | 90.0 % | 91.5 % |
| alternative, competing | 19 760 | 82.2 % | 88.0 % |

* **Contigs match reads except for minor isoforms.** 41 % of the 131 580
  contig junctions are annotated, 22 % of the 249 552 read junctions.
* **The loss follows relative usage.** Competing introns with under 5 % of
  the reads of their competitor: contigs find 71 %; dominant ones: 97 %.

### Unitigs keep the branches the contigs lose (10 runs)

Logan headers carry BCALM-style links (`L:+:89796:-`) in both files, but
the contigs have almost none left (SRR10620183: 92 % of contigs without a
link, against 14 % of unitigs), so the minor branch of a bubble is gone.
For 10 accepted runs (different BioProjects, spread over the ranking) the
junction 31-mers of every reference intron (≥ 10 nt on each side) were
looked up in the unitig and contig files:

| class | contigs | unitigs |
|---|---|---|
| single-isoform gene | 82.2 % | 88.0 % |
| constitutive | 94.4 % | 96.9 % |
| alternative, non-competing | 70.8 % | 82.9 % |
| alternative, competing | 58.0 % | 83.3 % |

* **Contig construction, not coverage, removes the alternative branch:**
  +2.5 points for constitutive introns, +25 for competing ones. Per run,
  unitigs contain 1.6–4.3× more competing junctions than contigs.
* 73 % of the 1 805 competing introns that A8's reads find and the
  375-run contig alignments miss have their junction k-mers in the
  unitigs of these 10 runs.

### Branch-path prototype (job 8232001, node236)

Per run: oriented unitig graph from the `L:` links; every node with ≥ 2
successors is a branch point; dead-end successors < 100 bp (error tips) and
short successors that rejoin at the same node with lengths within 1 bp
(SNP / 1-bp indel bubbles) are dropped; each remaining branch becomes a
path of ≤ 150 bp before and ≤ 150 bp after the branch point, extended
through the highest-abundance neighbours. Paths and contigs are aligned
with `minimap2 -ax splice --secondary=no -G 20000`; introns need ≥ 20 bp
aligned on both sides.

| class | contigs | + branches, canonical, ≥ 2 runs | + branches, canonical | k-mer ceiling |
|---|---|---|---|---|
| single-isoform gene | 82.7 % | 83.2 % | 84.2 % | 88.0 % |
| constitutive | 95.4 % | 95.6 % | 96.0 % | 96.9 % |
| alternative, non-competing | 67.2 % | 72.8 % | 77.9 % | 82.9 % |
| alternative, competing | 50.5 % | 65.8 % | 76.7 % | 83.3 % |

* **Branch paths recover 92 % of what the unitigs contain** for
  competing introns (76.7 of 83.3 %).
* **Unfiltered they are noisy.** 55 245 introns come only from branch
  paths; 81 % have a canonical splice site and 14 % of those are annotated.
  With ≥ 2 supporting runs 8 651 remain, 40 % annotated; 57 % of the novel
  ones are also in A8's reads (novel contig introns: 49 %). The 2-run rule
  keeps +15 points of competing introns at contig-like precision; with
  more runs it should keep more.
* **Cost.** Unitig files are 10–20× the contig files (these 10 runs:
  1.1 GB, one other accepted run 3.7 GB), so unitigs are worth fetching
  for a few dozen top runs, not for all accepted runs. Download on brain
  ~10 MB/s. The Python extraction is the bottleneck: 18 s to 25 min per
  run, 0.8–3.8 M paths (up to 1 GB) for the larger runs; the filters
  still pass most error branches and need tightening (e.g. require
  unequal branch lengths, drop low-abundance branches) before this goes
  into `varus logan`. Alignment of the paths: 7–235 s per run.

### *Takifugu* end-to-end test and decision (2026-09-25)

Full chain (Logan contigs of 300 accepted runs, optionally plus filtered
branch paths from the unitigs of 30 runs → minimap2 splice → StringTie 3
`-L` → DRUSILLA vertebrates, `--lorf-class --drop-unstranded`) scored as
in the GCB poster (gffcompare `--strict-match -e 3`, CDS level) on
*Takifugu rubripes*:

| Evidence | gene F1 (S / P) | tx F1 | + Tiberius (tib_chp) |
|---|---|---|---|
| VARUS reads (poster baseline) | 71.59 (63.0 / 82.9) | 50.87 | 80.05 |
| Logan contigs | 62.59 (57.9 / 68.1) | 39.44 | 76.54 |
| Logan contigs + branches | 63.62 (59.1 / 68.9) | 41.40 | 76.80 |
| Logan contigs + branches, StringTie `--min-cov 2` | 63.75 (59.1 / 69.2) | 41.51 | 76.79 |

The gap is precision. Coverage/TPM filters before or after DRUSILLA do
not close it; gffcompare class "j" (partial junction match) roughly
doubles. Contigs appear to flatten isoform abundance, so the wrong
isoforms are picked. The pipeline is also not faster than VARUS
sampling once branch extraction is included.

**Decision: not pursued further.** Logan stays a prescreen for run
selection; transcript evidence for annotation keeps coming from sampled
reads. Scripts and outputs: brain
`~/varus_bench/logan_branches/takifugu`.

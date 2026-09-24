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
| A7_defaults | (v2 defaults, 2026-09-24: `--parallel-downloads 6`) | the shipped configuration without Logan |
| A8_logan_defaults | A7 + `varus logan` (16 download connections, `-K 20M`, 25-run chunks, streamed scan) + `--logan-dir` | the shipped configuration with Logan, before the sparse estimator |
| A9_defaults_sparse | A7 + sparse estimator | the shipped configuration without Logan |
| A10_logan_sparse | A8 + sparse estimator; Logan back to 8 connections, no `-K` | the shipped configuration with Logan |

A3/A3t2 used `--prefetch`, which was retired on 2026-09-24; they are kept
for the record and were removed from the script. A1/A2/A4–A6 now pass
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

* Coelastrella tenuitheca GCA_051903525.1 — 9 runs, no rejections (tests
  the in-loop speed-ups in isolation; Logan can prune nothing here).
* Chlorella sorokiniana GCA_025917655.1 — 391 runs, 35 % of batches rejected
  in production (tests the Logan gate).
* Drosophila melanogaster GCF_000001215.4 — 115 457 runs, RefSeq annotation
  (tests the many-run case and intron Sn/Sp against a real annotation, as in
  the VARUS paper).

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

Chlorella sorokiniana (391 runs):

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

Coelastrella tenuitheca (9 runs):

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
  Tenuitheca), the per-batch BAM scan fell from 1.89 h to 10 min and the
  HISAT2 time from 1.74 h to 47 min. At 48 threads A1 equals A0, because
  alignment and scan are already small next to the downloads.
* **Downloads dominate everything else.** A0 spends 3.7–3.9 h of its
  4.5 h in `fastq-dump`. Three parallel downloads (A2) give 2.8–3.5×.
  Remote range dumps also get slower deep into a large run: Tenuitheca
  A2's last 100 batches took 30 s median.
* **Prefetch pays off when few runs are sampled many times, and costs time
  when many runs are sampled once.**
  * Tenuitheca: local `.sra` dumps take 1.3 s, so prefetch adds another
    2.5× on top of A2 (A3: 7.0× vs A0, 13.2× vs production).
  * Sorokiniana: A3 is slower than A2 (1.84 h vs 1.27 h). It prefetched
    93 GB, 18 GB of it for runs whose first batch then failed the
    quality gate, and the prefetches compete with the remote dumps for
    bandwidth.
* **Logan did not speed up Sorokiniana, and the benchmark shows why.**
  Only 17 of the 391 runs are usable (SRR37043103–123, 79–86 % unique
  HISAT2 alignments). All 17 are newer than the last Logan rebuild, so
  they are "absent". Nearly all other runs map below 5 % and are rejected
  after one batch each, which is 37 % of all batches.
  * Logan accepted 321 of these runs. Their contigs still cover the
    genome under `minimap2 -x splice` (2 300–4 600 tiles, yield 10–60 %).
  * Contig alignment identity shows they come from other species: the
    median `de` divergence of the contigs is **0.13–0.17**, with 0 % of
    aligned bases within 2 %. Tenuitheca's own runs score 0.0000 (100 %
    of bases within 1 %).
  * The breadth gate therefore ranked other-species runs first, and they
    took the bootstrap picks.
  * Fix: `varus logan --max-divergence` (default 0.05, see below).
* **Tenuitheca has nothing for Logan to remove.** All 9 runs are good, and
  the pre-screen takes 50 s. A4 and A3 are equal within noise (0.70 h vs
  0.65 h).
* **Half the batches cost ~10 % of the score.** A5 (500 batches) reaches
  89–90 % of A0's S on both species, below the 95 % rule. For Sorokiniana
  it spent 187 of its 500 batches on rejected runs.
* **`--profit-condition` never fired.** A6 ran all 1000 batches on both
  species; its differences from A4 are noise.
* **Network speed changes over the day.** Sorokiniana A0 downloads took
  9.1 s median against about 5.4 s in production, so A0 comes out 4 %
  slower than production despite 48 threads. Compare arms with A0
  (same day, same infrastructure) rather than with production.
* **greif14 (no SLURM, office uplink):** Sorokiniana A4 took 15.3 h
  (Logan 32 min). The run is network-bound on that machine.

### Rerun with the divergence gate and with 6 downloads (2026-09-24)

A4 and A5 on Sorokiniana were rerun with `varus logan --max-divergence
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
* **With the gate, A5 meets the old rules on Sorokiniana:** 95.2 % of A0's
  score in 16 % of its wall time. It recovers 92.5 % of A0's
  well-supported introns, against 94.8 % for a replicate. On Tenuitheca,
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

### Drosophila melanogaster (2026-09-24)

115 457 RNA-seq runs, 48 threads, all arms on the same day. The Logan arms
used the divergence gate. A4/A5 were submitted with `--prefetch` before it
was ruled out; their `varus run` time includes the prefetches. Intron
Sn/Sp use the VARUS paper's definition (Stanke et al. 2019): predicted =
distinct introns from the spliced alignments, 32 bp–350 kb (the `bam2hints`
window); reference = the 47 911 coding introns of the RefSeq annotation.
The paper reports Sn 0.935 / Sp 0.359 for Drosophila from 758 runs.

| arm | batches | rejected | wall | of which Logan | S / A0 | tiles ≥ 10 | runs sampled | Sn / Sp (all) | Sn / Sp (≥ 2 reads) | vs A0 |
|---|---|---|---|---|---|---|---|---|---|---|
| A0_baseline_flags | 1000 | 83 | 4.33 h | – | 100.0 % | 28979 | 265 | 0.950 / 0.199 | 0.933 / 0.293 | 1.0× |
| A2_parallel3 | 1000 | 86 | 1.58 h | – | 83.6 % | 26963 | 310 | 0.957 / 0.223 | 0.943 / 0.342 | 2.7× |
| A3_parallel3_prefetch | 1000 | – | loop 1.50 h, then hung | – | – | – | – | – | – | – |
| A4_logan_1000 (+prefetch) | 1000 | **27** | 3.94 h | 23 min | 92.8 % | 28123 | 198 | 0.957 / 0.184 | 0.936 / 0.301 | 1.0× |
| A5_logan_500 (+prefetch) | 500 | **5** | 2.65 h | 23 min | 83.6 % | 27647 | 118 | 0.930 / 0.298 | 0.896 / 0.446 | 1.4× |
| A7_defaults (6 downloads) | 1000 | 76 | **0.97 h** | – | 86.1 % | 27547 | 233 | 0.951 / 0.245 | 0.934 / **0.378** | **4.5×** |
| A8_logan_defaults (6 downloads, no prefetch) | 1000 | **18** | 2.36 h | 15.7 min | 92.7 % | 28091 | 165 | 0.954 / 0.186 | 0.932 / 0.305 | 1.8× |

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
* **A2's score is 16 % below A0 here, against 0.4 % on Sorokiniana.** A0
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

### Measured vs expected

| | expected (plan) | measured |
|---|---|---|
| P1 alone, 2 threads | 8.6 → 6.6 h | scan 1.89 h → 10 min, align 1.74 h → 47 min (Tenuitheca) |
| P1 + P2 + P3, 2 threads | 3–4× (8.6 → ~2.2 h) | Tenuitheca 7.1× (1.21 h); Sorokiniana 2.2× (1.91 h) |
| + real thread count | ~6× | Tenuitheca 13.2× (A3); Sorokiniana 3.4× (A2) |
| + Logan | −35 % batches where foreign runs exist | rejections 371 → 246; good runs absent from Logan, other species accepted |
| + half the batches | ÷2 if the score holds | ÷1.8, but S 89–90 % (fails the 95 % rule) |
| 6 downloads as default (Drosophila, 48 threads) | faster than 3 | 4.5× over A0, 1.6× over 3 downloads; S 86 %, Sn/Sp unchanged or better |
| Logan: SAM piped into the scanner, 16 connections, 25-run chunks | stage −1.5 to −3 min | stage +2.3 min (13.4 → 15.7 min); throughput unchanged at ~8 MB/s |
| Logan + 6 downloads (Drosophila) | −rejected batches, faster loop | rejected 76 → 18, S 86 → 93 %, but 2.4× slower than A7 (estimator 5.2 s/batch) |

### Decisions

* **Do not lower `--max-batches` in BRAKER4.** A5 fails the score rule on
  both species.
* **BRAKER4 wrapper defaults:**
  * `--parallel-downloads` defaults to 6 since 2026-09-24 (Sorokiniana:
    1.46× over 3 for 20 more rejected batches and 0.5 % of the score), so
    the wrapper passes nothing; K=1 still reproduces the v1 pick sequence;
  * pass the real thread count to the wrapper;
  * do **not** set `--prefetch` (decided 2026-09-24). It wins only for
    species with a handful of runs (Tenuitheca). With many runs it fills
    node-local disk with `.sra` files of runs that are then rejected
    (Sorokiniana 93 GB; watchdog stop on node385), competes with the range
    dumps for bandwidth, and prefetches queued near the end keep
    downloading after the last batch: Drosophila A3 finished its 1000
    batches in 1.5 h (same as A2) and then hung for 3 h in 58 queued
    prefetches until it was cancelled. The flag stays for experiments.
* **Use Logan with the divergence gate at 1000 batches.** On Sorokiniana
  it removed 97 % of the wasted batches (4.1× faster than A0, score
  105 %). On Tenuitheca it costs 50 s. Its value is removing runs from the
  wrong species, which the breadth gate alone could not do. On Drosophila
  it cut rejected batches from 83 to 27 and the intron sensitivity against
  the RefSeq annotation matches the VARUS paper. The stage costs 13–16 min
  there (23 min before `--scan-workers`). With the shipped defaults,
  though, the Drosophila run is 2.4× slower with Logan than without (A8
  2.36 h against A7 0.97 h) because the estimator limits the loop. The
  gain is score (93 % against 86 % of A0), not time. On species without
  foreign runs, Logan is a speed-up only once the estimator is faster.
* **Keep 1000 batches.** 500 batches meets the score rule only when Logan
  removes waste (Sorokiniana 95.2 %), not in general (Tenuitheca 89 %).
* **Replace the intron rule.** "≥ 90 % intron Jaccard" becomes "intron
  Jaccard at ≥ 5 reads no lower than the replicate floor (production vs
  A0)".

The raw tables come from `benchmark_varus.py report --outdir runs/<sp>
--baseline-log data/<sp>/baseline_varus.log`. The report reads
`baseline_Coverage.csv` and `baseline_introns.gff` next to the log.

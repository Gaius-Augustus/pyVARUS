# Changelog

All notable changes to VARUS are documented in this file.

## [2.0.0a0] -- unreleased

Full Python rewrite of the C++/Perl implementation. The online sampling
algorithm (Stanke et al., 2019,
[DOI:10.1186/s12859-019-3182-x](https://doi.org/10.1186/s12859-019-3182-x))
is unchanged.

### Changed vs. v1 (upstream [Gaius-Augustus/VARUS](https://github.com/Gaius-Augustus/VARUS))

| | v1 | v2 |
|---|---|---|
| Language | C++ + Perl + Bash | Python 3.9+ |
| Aligner | STAR or HISAT2 | HISAT2 (short reads), minimap2 (long reads, `--longreads`) |
| Read download | `fastq-dump --fasta` | `fastq-dump` (per-batch spot ranges) |
| Alignment intermediate | SAM -> samtools sort -> BAM | piped -> coordinate-sorted BAM directly |
| Intron extraction | `bam2hints` (AUGUSTUS) | `pysam` reimplementation |
| Strand assignment | `filterIntronsFindStrand.pl` | `pyfaidx` reimplementation |
| Final merge | hierarchical bash scripts | `samtools merge` |
| Pipeline driver | `runVARUS.pl` + `VARUSparameters.txt` | Nextflow + `varus` Python CLI |
| Per-iteration coverage dump | always (~8 GB for 1000 batches) | off by default, `--coverage-trace N` |
| Per-batch FASTA kept gzipped | yes | deleted by default, `--keep-batches` to retain |
| User-facing parameters | ~25 in a parameters file | ~10 CLI flags + `--advanced KEY=VALUE` |

### Added (speed-ups, 2026-09)

- `varus logan`: optional pre-screen that aligns each candidate run's Logan
  contigs (public S3) to the genome, rejects foreign/empty runs by tile
  breadth and other species by contig divergence (`--max-divergence`), ranks the rest by the VARUS score, writes `Runlist.logan.tsv`,
  `logan/LoganRanking.tsv`, `logan/logan_introns.gff` and a seed splice DB.
  `varus run --logan-dir` consumes it (run filter, estimator prior,
  splice-DB seed). Requires minimap2 and the `[logan]` extra (`zstandard`).
- `varus logan`: minimap2's SAM is piped straight into `--scan-workers`
  (default 2) scanner processes, so a chunk is scanned while it aligns and
  no chunk BAM is sorted or written unless `--logan-bam` is set. Before,
  the single-threaded scan of each chunk BAM ran after its alignment and
  was as slow as the alignment (Drosophila: 12 min of a 23 min stage);
  scanning from disk in worker processes brought the stage to 13 min; the
  streamed scan removes the sort and the disk round trip. `--chunk-runs`
  defaults to 25 (each minimap2 call reloads the index). 16 download
  connections and `-K 20M` were tried and dropped: the S3 rate stayed at
  ~8 MB/s in total and minimap2 got slower (see `docs/benchmark_logan.md`).
- `varus run --parallel-downloads K` (default 6): K concurrent batch
  downloads with lazy-greedy picks that account for in-flight batches'
  expected gains. `fastq-dump` on a remote spot range is latency-bound, so
  K=3 gave 2.7–3.5× and K=6 4.5–5× over serial downloads on the benchmark
  genomes (Drosophila: 4.33 h → 58 min). K=1 reproduces the v1 pick
  sequence exactly.
- `varus run` aligns the next downloaded batch in a background thread while
  the main thread scans and scores the current one (on whenever downloads
  are pipelined; `--no-align-ahead` turns it off). Picks are unchanged; the
  aligner may use a splice DB that lacks the current batch's junctions. The
  splice DB is now replaced atomically. `BatchTimings.tsv` gains `t_wait`
  (main thread blocked on download or alignment); `t_align` now overlaps
  other phases, so the phase columns no longer add up to `t_batch`.
- `varus run --merge-batches N` (default 10; 1 = off; not used with
  `--parallel-downloads 1`, which keeps reproducing the v1 pick sequence):
  repeated picks of a run come as one download of a contiguous spot range
  of up to N batches, aligned and scanned once. An extra batch is claimed
  only while greedy selection would pick the same run again with the
  claimed batches counted as observed, and only after the run's first
  batch passed the quality gate. `--max-batches` still counts 50 k-spot
  batches. Each `fastq-dump` call has a fixed cost of 5–27 s (500 k spots
  of SRR36274151: 34 s against 27 s for 50 k), which dominated Logan runs
  that exploit one such run. Drosophila, 1000 batches: loop −37 % without
  Logan (29 → 18 min, S unchanged on average), −63 % with Logan and
  λ = 3 (55 → 20 min; 34 min including the Logan stage, S 101 % of A0).
  `BatchTimings.tsv` gains `n_batches`.
- `varus run --scan-workers N` (default 4 at 48 threads, see below; 0/1 = off): merged batches are
  indexed and scanned by genome region in N processes; names split across
  regions or with more than one hit are resolved centrally, so the counts
  are identical to the one-pass scan (checked on three real Drosophila
  batches). Merged-batch scan 3.1–6.8 → 1.8–3.5 s; on brain the loop is
  then download bound, so its wall time did not change.
- `varus logan --align-groups N` (default 3 at 48 threads, see below): each chunk is aligned by N
  minimap2 processes (threads split evenly), each piped into its own
  scanner, because one scanner throttled minimap2. Rankings are
  byte-identical; the Drosophila Logan stage takes 11.5 instead of
  14.3 min. Genomes over 1 Gb use one group (one index copy per process).
- Estimator defaults `lambda=3 pseudo-count=0.1` (v1: 10 and 1). Same or
  higher score in every benchmark, 3 seeds each at 1000 batches: Tenuitheca
  +0.3–0.6 %, Sorokiniana without Logan +1 % (99.7 against 98.7 % of A0),
  with Logan equal (106.2 %); Drosophila 76.6 → 91.7 % without Logan and
  93.2 → 102 % with it, because the loop stays on good runs instead of
  returning to the pooled profile. `--advanced lambda=10 pseudo-count=1`
  restores v1.
- Large runlists (mouse: 2 M runs, 860 M batch indices): each run's
  shuffled batch order is built on its first pick instead of at load time
  (one seed per run is drawn at load), which avoids ~35 GB and minutes of
  shuffling. A given `--seed` therefore samples different spot ranges within
  a run than earlier 2.0 alphas (same distribution). From `fresh_pool_min` (1000) never-downloaded runs without a
  Logan prior on, these interchangeable runs are kept as one pool: one
  member stands in for all of them in the estimator, the profit and the
  lazy greedy, and the avg_len-weighted tie-break draws over all of them in
  runlist order exactly as before. Picks are identical with and without the
  pool (tests: serial loop, lazy greedy with in-flight extras, parallel
  refill). The periodic `RunStatistics.csv` omits pool members (all-zero
  rows); the final one lists every run.
- Thread budget: stages that run at the same time now share `--threads`
  instead of each taking all of it. In `varus run` the aligner gets
  `--threads` minus one core for the main thread (align-ahead), one for the
  `fastq-dump` processes and the rolling merge's `-@` (at most 4) while a
  merge runs; `samtools sort -@` per batch is at most 4 (was threads − 1).
  Before, a 48-thread run started HISAT2 `-p 48`, `samtools sort -@ 47` and
  a rolling merge `-@ 24` beside the scan and six downloads. In
  `varus logan` minimap2 gets `--threads` minus the scanner processes and
  one core for the main and download threads. Reservations are capped at a
  quarter of `--threads`. Index builds and final merges run alone and keep
  all threads.
- `--scan-workers` (`varus run`) and `--align-groups` (`varus logan`) default
  to values derived from `--threads`: `threads // 8` scan workers, at most
  4, none below 16 threads; one minimap2 group per ~15 minimap2 threads,
  1–4. 48 threads give the benchmarked 4 and 3. Before, 4 scanner processes
  and 3 minimap2 processes (three copies of the index) also started on a
  4-core machine. Explicit values still win. README: "Threads and machine
  size"; figure `docs/figures/speedups.svg`.
- `varus run --prefetch`: `prefetch` a run's `.sra` after its second pick
  and range-dump locally (per-run and total disk caps, LRU eviction).
  Experimental and off by default; not recommended (fills local disk,
  queued prefetches outlive the batch loop; see `docs/benchmark_logan.md`).
- Rolling background merge of batch BAMs (`--merge-every`, default 100);
  the final merge only joins the parts.
- `BatchTimings.tsv` with per-batch phase timings; `TIMING` log lines.
- Exit status 3 when no batch passed the quality gate (previously exit 0
  with no `VARUS.bam`, which crashed downstream wrappers).
- `Runlist.tsv` gains a `bioproject` column (optional when reading).

### Changed (speed-ups, 2026-09)

- Estimator and profit work on sparse per-run counts against a tile index
  that only grows. Before, every batch rebuilt a dense array per run with
  data from a Python dict over all tiles, plus a `run.p` dict nobody read:
  with Logan's 375 priors on Drosophila (29 k tiles) that was 5.2 s per
  batch, 1.4 h of a 2.1 h run; now 0.2 s. `run.p` is no longer filled.
- Each batch BAM is scanned once for UMRs, spliced reads and introns
  (`tiles.scan_batch_bam`); before, UMR counting and intron extraction each
  read the whole BAM (0.40 → 0.28 s per 100 k records, same counts).
- Splice-site DB is maintained incrementally: introns are stranded once
  (`StrandAssigner` cache) and the DB is rewritten only when new junctions
  appear. Previously every batch re-stranded all cumulative introns, which
  grew from 2 s to 10 s per batch over a 1000-batch run.
- HISAT2 runs with `--mm --no-unal` and per-batch BAMs use compression
  level 1 (`--no-hisat2-mm`, `--keep-unaligned` restore the old behaviour).
  `VARUS.bam` therefore no longer contains unaligned reads.
- Rejected batches and failed downloads no longer leave BAM/log files or
  empty directories behind.
- Estimator accepts per-run pseudo-observations (Logan prior); with none
  given it is unchanged.

### Added

- Long-read RNA-seq support (`--longreads`): minimap2 alignment with per-run
  preset selection from SRA platform metadata; BED12 splice-DB feedback;
  uniqueness % from BAM scan instead of HISAT2 log parsing.
- Standalone Nextflow pipeline in [`nextflow/`](nextflow/).
- Entrez retry-with-backoff on HTTP 429 / 5xx in `varus runlist`.

### Fixed (2026-09-25)

- `varus --version` (documented in CONTRIBUTING, was missing).
- The `species` argument of `varus run` is logged at start (it was parsed
  and ignored); the per-batch `TIMING` log line reports the score S.
- Nextflow `VARUS_LOGAN`: `varus logan` exit 3/4 aborted the process under
  `set -e` although the comment promised a fallback. Exit 3 now keeps the
  runs Logan could not screen (without the prior), exit 4 the full runlist;
  `nextflow/NO_LOGAN/` exists so `main.nf` can stage it when Logan is off.
- Dead code removed: `RunState.p` and its fallback loops, `total_profit`,
  `_pick_and_download_single`, unused imports and constants.

### Removed

- Per-iteration coverage dump is now opt-in (`--coverage-trace N`).
- Per-batch FASTA/BAM are deleted after counting (`--keep-batches` to retain).

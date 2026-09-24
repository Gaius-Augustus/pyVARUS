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
| Read download | `fastq-dump --fasta` | `fastq-dump` (per-batch ranges) + `fasterq-dump` (full runs) |
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

### Removed

- Per-iteration coverage dump is now opt-in (`--coverage-trace N`).
- Per-batch FASTA/BAM are deleted after counting (`--keep-batches` to retain).

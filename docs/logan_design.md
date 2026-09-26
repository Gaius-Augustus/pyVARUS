# Design: Logan-assisted pyVARUS

This document records why and how Logan was integrated into pyVARUS. It
supersedes the earlier brief `logan-instructions/VARUS_LOGAN_IMPLEMENTATION.md`,
which proposed a stand-alone replacement package; the decision (2026-09-23)
was to *speed up* pyVARUS instead of replacing it.

## Where pyVARUS spends its time

Two production runs on the group cluster (BRAKER4 wrapper, 2 threads, 50 k
read pairs per batch) were profiled from their logs:

| | *Coelastrella tenuitheca* (9 runs) | *Chlorella sorokiniana* (391 runs) |
|---|---|---|
| Loop wall time | 8 h 15 min (1000 batches) | 4 h 09 min (~1000 batches) |
| Per batch | 29.7 s | 15.8 s |
| `fastq-dump` (remote spot range) | 16.6 s (56 %); 8.5 s early, ~28 s late | 8.5 s (54 %) |
| HISAT2 + samtools sort | 6.3 s (21 %); 3.6 → 10 s as the splice DB grew | 3.9 s (25 %) |
| BAM scan + introns + strand DB + estimator | 6.8 s (23 %); 2.2 → 9.9 s | 3.4 s (21 %) |
| Batches rejected by the quality gate | 0 | 343 (≈35 %) |
| Final `samtools merge` | 22 min | 8 min |

The loop was serial and latency-bound, two per-batch costs grew with the
number of introns seen, and every foreign run cost one full batch before
the gate rejected it.

## What Logan offers

[Logan](https://github.com/IndexThePlanet/Logan) holds unitigs and contigs
(k = 31, k-mers seen once discarded) for nearly every public SRA run on
`s3://logan-pub`, public and credential-free. It is rebuilt in full at
discrete time points (last ≈ Aug 2026), not incrementally, so availability is
decided per run by an HTTP `HEAD`.

Measured on *S. pombe* (2026-09-23):

* contigs are 2.7–5 MB (zstd) per RNA-seq run; unitigs 16–420 MB;
* `minimap2 -ax splice --secondary=no` on one run's contigs: 3.6 s;
  4 250 introns, precision 0.93 / sensitivity 0.75 against PomBase,
  2 436 of 2 514 5-kb tiles covered; unitigs gave sensitivity 0.56 at 2× cost;
* 52 runs: 11 s of downloads, 118 s of alignment; foreign or empty runs cover
  < 10 % of the tiles of a real run; a greedy pick of 5 runs already covers
  2 499 tiles and 88 % of reference introns;
* 83 % of a random sample of *S. pombe* runs had Logan contigs.

## Decisions (and what changed relative to the brief)

| Brief | Decision | Reason |
|---|---|---|
| New package `logan_varus` | `varus logan` stage + `varus run --logan-dir` inside pyVARUS | keep algorithm, outputs and BRAKER4 interface |
| Seeds from Vipsania / miniprot + junction k-mer counting | none; align contigs to the genome | minimap2 is already a dependency; gives tiles *and* introns in seconds; no gene predictions, no compiled k-mer tool |
| Rank on unitigs | rank on contigs | 100× smaller, higher intron sensitivity in the test |
| "Logan ends 2023" | freeze-based availability, per-run HEAD, read fallback | verified 2025 runs present |
| Greedy set cover over junctions | greedy on the VARUS score Σ log(1 + c_j) | same objective as the online loop |
| Gate on mapped fraction | gate on tile breadth | an rRNA-heavy run had 7 % mapped contigs but 2 410 tiles |
| Contig BAM merged into the evidence BAM | `LOGAN.bam` optional and separate | BRAKER4 consumes read alignments only |
| Rank by breadth only | breadth × yield, yield = mapped fraction of abundance × length contig mass | breadth alone ranked a mixed sample first (2.6 % of contigs mapped, 5.5 % unique reads); the contig-*count* mapped fraction does not predict read yield (a 5.7 %-mapped run gave 85 % unique reads), abundance × length does separate good (≥ 35 %) from bad (≤ 19 %) runs on 9 calibration runs |
| Estimator prior alone decides first picks | Logan bootstrap: first picks follow the ranking | at cold start the uniform shared prior out-scores every informative profile (Jensen), so never-processed runs were picked first |

## How the pieces fit

```
varus runlist  ──► Runlist.tsv (8 columns, incl. bioproject)
varus index    ──► HISAT2 index
varus logan    ──► HEAD ► contigs ► minimap2 (chunks) ► tiles + introns per run
                   ► breadth gate ► greedy ranking
                   ► Runlist.logan.tsv, logan/LoganRanking.tsv,
                     logan/logan_introns.gff, logan/logan.splice_sites,
                     logan/logan_tiles.tsv.gz, logan/logan_summary.json
varus run --logan-dir
               ──► drops rejected runs, seeds intronDB before batch 1,
                   first picks = ranked runs in rank order (bootstrap),
                   estimator prior  p̂_r ∝ c^r + ℓ^r + a + λ·T·p̄
                   (ℓ^r = contig tile profile scaled to --logan-prior-batches
                   batches; ℓ = 0 reproduces v1 exactly),
                   expected reads of a run × its yield (accepted) or × the
                   gate's acceptance rate (unprocessed runs: kept only if
                   the accepted ones hold < --max-batches batches, or with
                   --logan-keep-unprocessed)
```

In-loop speed-ups independent of Logan:

* `StrandAssigner` strands each junction once; the splice-site DB is
  rewritten only when new junctions appear.
* `--parallel-downloads K`: K downloads in flight; picks are lazy greedy
  against `total_obs + Σ expected(in-flight)`. With K = 1 the pick sequence
  is bit-identical to v1 (verified against the pre-change controller over
  5 seeds × 40 batches).
* Align-ahead: the next finished download is aligned in a one-thread
  executor while the main thread scans and scores the current batch; the
  batch keeps its in-flight accounting until it is applied, and the splice
  DB is replaced atomically.
* `--prefetch` (`prefetch` a run's `.sra` after its second pick, local
  `fastq-dump -N/-X` afterwards) was implemented, benchmarked and removed
  (2026-09-26): it fills local disk with runs that are then rejected and
  keeps downloading after the loop; merged batches address the same cost.
* HISAT2 `--mm --no-unal`, level-1 per-batch BAMs, rolling background merge.
* Exit status 3 when nothing passed the gate.

## Expected effect

From the measured profile (*C. tenuitheca*-like, 2 threads): incremental DB and
merge changes ≈ 8.6 h → 6.6 h; parallel downloads ≈ 2.1 h; prefetch ≈ 2.0 h
(since removed);
the real thread count from BRAKER4 ≈ 1.3 h; Logan removes the 35 % of
batches spent on foreign runs where they exist and makes the first picks
informed. Fewer batches for the same score is a policy decision validated by
`docs/benchmark_logan.md`.

## Known limits

* Logan contig introns favour well-expressed genes (k-mers seen once are
  dropped); read batches still add the low-coverage junctions.
* Runs released after the last Logan rebuild are "unprocessed" and, since
  2026-09-26, dropped by default (`--logan-keep-unprocessed` keeps them with
  the shared prior). They are kept automatically when the accepted runs
  hold fewer batches than `--max-batches` or when Logan accepted no run.
* Cross-strain runs can pass the breadth gate yet fail HISAT2's 5 % gate;
  the in-loop gate still catches them.
* `varus logan` caps candidates (default 500, round-robin over BioProjects)
  for species with tens of thousands of runs.

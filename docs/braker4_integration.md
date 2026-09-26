# Using the pyVARUS v2 speed-ups from BRAKER4

pyVARUS v2 keeps the `runlist → index → run` interface BRAKER4 already uses
(`scripts/run_varus_wrapper.sh`) and its outputs (`VARUS.bam`, `introns.gff`,
`Coverage.csv`, `RunStatistics.csv`). Nothing in BRAKER4 *must* change. The
diffs below switch on the speed-ups and fix the known "no usable RNA-seq"
crash. They are listed here rather than applied, because BRAKER4 lives in a
separate repository.

## 1. `docker/pyVARUS/Dockerfile`

No change needed. The current image installs everything pyVARUS v2 uses:

```dockerfile
RUN micromamba install -y -n base \
        -c conda-forge -c bioconda \
        python=3.11 \
        hisat2 \
        minimap2 \
        "samtools>=1.17" \
        sra-tools \
        pip \
        git && \
    micromamba clean -afy
RUN git clone https://github.com/Gaius-Augustus/pyVARUS.git /opt/pyvarus && \
    pip install --no-cache-dir -e "/opt/pyvarus[align]"
```

The image already has `minimap2`, which the Logan pre-screen needs; the
`zstandard` package that decompresses the contigs is a core dependency of
pyVARUS, so the `pip` line is unchanged. Bump the tag
(`docker/pyVARUS/build_push.sh`, `Snakefile` `varus_image`,
`config.ini.example`) e.g. to `katharinahoff/pyvarus:v2.0.0`.

## 2. `scripts/run_varus_wrapper.sh`

The Logan pre-screen is part of `varus run` (on by default since
2026-09-25): it screens the runs before the first download, writes
`<outdir>/logan/` and `Runlist.logan.tsv`, samples only the accepted runs
(the runs Logan could not screen are added when the accepted ones hold
fewer batches than `--max-batches`, or with `--logan-keep-unprocessed`),
and falls back by itself when Logan accepts no run (exit 3: the runs Logan
could not screen are sampled without a prior) or is unreachable (exit 4:
the loop runs without the pre-screen). Nothing has to be added between step 2
(index) and step 3 (run). Step 3 only needs to handle exit status 3 of
`varus run`:

```bash
set +e
varus run "$SPECIES_NAME" "$GENOME_ABS" \
    --runlist "$VARUS_DIR_ABS/Runlist.tsv" \
    --index   "$INDEX_DIR/hisatidx" \
    --outdir  "$VARUS_DIR_ABS" \
    --threads "$THREADS" \
    >> "$LOGFILE_ABS" 2>&1
rc=$?
set -e
if [ "$rc" = "3" ]; then
    # No batch passed the quality gate: no usable public RNA-seq for this
    # genome. Make it a data condition, not a crash.
    echo "[ERROR] pyVARUS found no usable RNA-seq for $SPECIES_NAME" >> "$LOGFILE_ABS"
    touch "$VARUS_DIR_ABS/NO_RNASEQ"
    exit 3
fi
[ "$rc" = "0" ] || exit "$rc"
```

To run the pre-screen with non-default options (`varus logan --help-all`),
call `varus logan ... --outdir "$VARUS_DIR_ABS"` before step 3; `varus run`
then reuses `$VARUS_DIR_ABS/logan/`. `--no-logan` skips the pre-screen.

Step 4: `VARUS.bam` is already coordinate-sorted (every batch BAM is sorted
and `samtools merge` preserves order), so the extra `samtools sort` pass can
be replaced by a header check:

```bash
if samtools view -H "$VARUS_DIR_ABS/VARUS.bam" | grep -q 'SO:coordinate'; then
    mv "$VARUS_DIR_ABS/VARUS.bam" "$OUTPUT_BAM_ABS"
else
    samtools sort -@ "$THREADS" -o "$OUTPUT_BAM_ABS" "$VARUS_DIR_ABS/VARUS.bam"
fi
samtools index -c -@ "$THREADS" "$OUTPUT_BAM_ABS"
```

## 3. `rules/preprocessing/run_varus.smk`

* Cleanup: also remove `logan/contigs`, `logan/bams` and `merged/` (all
  scratch). Keep `logan/LoganRanking.tsv`, `logan/logan_summary.json`
  and `BatchTimings.tsv` next to the other diagnostics.
* Threads: the production logs (all 8 `logs/*/varus/varus.log` of the
  chlorophyte run, 2026-09-12) show `[INFO] Threads: 2` although
  `config.ini` sets `cpus_per_task = 48`. Cause (found 2026-09-25): the
  project scripts `bulk_algae_braker4.sh` and `chloro_annotation.sh` start
  Snakemake with `--cores 1 --executor slurm`. Snakemake 9.24 caps every
  rule's `threads:` at `--cores`, so the SLURM plugin submitted all 48 jobs
  of that run with `--cpus-per-task=1` (sacct: AllocCPUS 1). Brain nodes
  have 2 hardware threads per core, which is where the 2 comes from. This
  applies to every rule, not only VARUS. Fix: `--cores 48` (matching
  `cpus_per_task`), as the BRAKER4 README already shows. The rule itself is
  fine; with 2 threads HISAT2 alone costs 4–10 s per batch.
* Optional: expose `--max-batches` in `config.ini` so it can be lowered once
  the benchmark (`docs/benchmark_logan.md`) confirms the score/intron parity
  of Logan-informed runs at 500 batches.
* Optional per-sample fallback: when `varus/NO_RNASEQ` exists, route the
  sample to protein-only (EP) mode instead of failing the DAG.

## 4. What the outputs mean now

| File | Change |
|---|---|
| `VARUS.bam` | unchanged format; unaligned reads are no longer included (`--no-unal`) |
| `introns.gff` | unchanged (read-derived). Logan contig introns are in `logan/logan_introns.gff` and only merged in with `--logan-merge-introns` |
| `BatchTimings.tsv` | new: per-batch phase timings (download / align / scan / DB / estimator) |
| `Runlist.logan.tsv` | new: the runlist after the Logan gate, ranked |
| `logan/LoganRanking.tsv` | new: per-run contig alignment breadth, gate status, greedy rank |
| exit status 3 | new: no batch passed the quality gate; `VARUS.bam` absent |
| exit status 4 (`varus logan`) | new: Logan S3 unreachable |

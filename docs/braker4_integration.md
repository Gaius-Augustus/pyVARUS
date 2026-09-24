# Using the pyVARUS v2 speed-ups from BRAKER4

pyVARUS v2 keeps the `runlist → index → run` interface BRAKER4 already uses
(`scripts/run_varus_wrapper.sh`) and its outputs (`VARUS.bam`, `introns.gff`,
`Coverage.csv`, `RunStatistics.csv`). Nothing in BRAKER4 *must* change. The
diffs below switch on the speed-ups and fix the known "no usable RNA-seq"
crash. They are listed here rather than applied, because BRAKER4 lives in a
separate repository.

## 1. `docker/pyVARUS/Dockerfile`

```diff
 RUN micromamba install -y -n base \
         -c conda-forge -c bioconda \
         python=3.11 \
         hisat2 \
         minimap2 \
         "samtools>=1.17" \
         sra-tools \
+        zstd \
         pip \
         git && \
     micromamba clean -afy
@@
-RUN git clone https://github.com/Gaius-Augustus/pyVARUS.git /opt/pyvarus && \
-    pip install --no-cache-dir -e "/opt/pyvarus[align]"
+RUN git clone https://github.com/Gaius-Augustus/pyVARUS.git /opt/pyvarus && \
+    pip install --no-cache-dir -e "/opt/pyvarus[align,logan]"
```

Bump the tag (`docker/pyVARUS/build_push.sh`, `Snakefile` `varus_image`,
`config.ini.example`) e.g. to `katharinahoff/pyvarus:v2.0.0`. The `logan`
extra only adds the pure-Python `zstandard` package; `zstd` (CLI) is a
fallback decompressor.

## 2. `scripts/run_varus_wrapper.sh`

Between step 2 (index) and step 3 (run):

```bash
# Step 2b: Logan pre-screen (optional; skipped when Logan is unreachable)
echo "[INFO] Logan pre-screen..." >> "$LOGFILE_ABS"
LOGAN_ARGS=""
set +e
varus logan "$GENOME_ABS" \
    --runlist "$VARUS_DIR_ABS/Runlist.tsv" \
    --outdir  "$VARUS_DIR_ABS" \
    --threads "$THREADS" \
    >> "$LOGFILE_ABS" 2>&1
rc=$?
set -e
case "$rc" in
  0) LOGAN_ARGS="--logan-dir $VARUS_DIR_ABS/logan"; RUNLIST="$VARUS_DIR_ABS/Runlist.logan.tsv" ;;
  3) # No run Logan could screen passed its gate (e.g. every run in Logan is
     # another species). Runlist.logan.tsv still lists the runs Logan could
     # not screen (too new for the last Logan rebuild), which are often the
     # only usable ones; use the full runlist only if it is empty.
     if grep -qv '^#' "$VARUS_DIR_ABS/Runlist.logan.tsv" 2>/dev/null; then
       echo "[WARN] Logan accepted no run; using the runs Logan could not screen" >> "$LOGFILE_ABS"
       RUNLIST="$VARUS_DIR_ABS/Runlist.logan.tsv"
     else
       echo "[WARN] Logan accepted no run; falling back to the full runlist" >> "$LOGFILE_ABS"
       RUNLIST="$VARUS_DIR_ABS/Runlist.tsv"
     fi ;;
  4) echo "[WARN] Logan unreachable; falling back to the full runlist" >> "$LOGFILE_ABS"; RUNLIST="$VARUS_DIR_ABS/Runlist.tsv" ;;
  *) echo "[ERROR] varus logan failed (rc=$rc)" >> "$LOGFILE_ABS"; exit "$rc" ;;
esac
```

Step 3 becomes:

```bash
set +e
varus run "$SPECIES_NAME" "$GENOME_ABS" \
    --runlist "$RUNLIST" \
    --index   "$INDEX_DIR/hisatidx" \
    --outdir  "$VARUS_DIR_ABS" \
    --threads "$THREADS" \
    --parallel-downloads 3 --prefetch \
    $LOGAN_ARGS \
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

* Cleanup: also remove `logan/contigs`, `logan/bams`, `sra/` and `merged/`
  (all scratch). Keep `logan/LoganRanking.tsv`, `logan/logan_summary.json`
  and `BatchTimings.tsv` next to the other diagnostics.
* Threads: the two production logs analysed on 2026-09-23 show
  `[INFO] Threads: 2` although `config.ini` sets `cpus_per_task = 48`. The
  rule's `threads:` must receive the real allocation; with 2 threads HISAT2
  alone costs 4–10 s per batch.
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

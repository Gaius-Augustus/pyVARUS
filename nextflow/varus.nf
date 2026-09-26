// VARUS v2 Nextflow module — four reusable processes that wrap the Python
// CLI. Designed to be `include`d from a parent workflow:
//
//     include { VARUS_RUNLIST; VARUS_INDEX; VARUS_LOGAN; VARUS_RUN } from '/path/to/pyVARUS/nextflow/varus.nf'
//
// Each process expects the `varus` CLI on $PATH (`pip install -e .[align]`)
// plus `hisat2`, `hisat2-build`, `samtools`, `fastq-dump` and `minimap2`.
// With `--longreads` set, `minimap2` replaces hisat2/hisat2-build.
//
// The processes feed each other:
//
//     VARUS_RUNLIST -> Runlist.tsv (NCBI Entrez query for the species)
//     VARUS_INDEX   -> HISAT2 (or minimap2 with --longreads) index of the genome
//     VARUS_LOGAN   -> pre-screen from Logan contigs (on by default; --varus_logan false)
//     VARUS_RUN     -> online loop: download SRA batches, align, score tiles
//
// Inputs are passed as a single tuple beginning with `species` and `genome`;
// callers may extend the tuple with arbitrary trailing fields, which are
// forwarded unchanged to make plumbing into larger pipelines easy.

// Memory and time scale with the genome FASTA (bytes ~ bases). Measured on
// wheat (14.8 GB FASTA, 2026-09-26): hisat2-build 121.5 GB peak (8.2 x the
// FASTA), minimap2 index build 67 GB per 8 Gbp part and 88.5 GB in one part,
// HISAT2 alignment 26 GB, minimap2 alignment ~30 GB per 8 Gbp part, and the
// `varus logan` main process ~25 GB for 500 runs. A task killed for memory
// (exit 137-140) is retried with twice the memory and time.
def genomeGb(genome) { genome.size() / 1e9 }
def memGb(double gb, int attempt) { "${(long) Math.ceil(gb * attempt)} GB" }
def hours(double h, int attempt) { "${(long) Math.ceil(h * attempt)}h" }

process VARUS_RUNLIST {
    tag { species }
    publishDir { "${params.outdir}/${species.replaceAll(' ', '_')}/varus" }, mode: 'copy', overwrite: true
    cpus 1

    input:
        tuple val(species), path(genome), val(extra)

    output:
        tuple val(species), path(genome), path("Runlist.tsv"), val(extra)

    script:
    def maxRuns  = params.varus_max_runs ?: 0
    def email    = params.ncbi_email ?: ''
    def apiKey   = params.ncbi_api_key ?: ''
    def emailArg = email  ? "--email ${email}"     : ''
    def keyArg   = apiKey ? "--api-key ${apiKey}"  : ''
    def longArg  = params.longreads ? '--longreads' : ''
    """
    set -euo pipefail
    varus runlist '${species}' \\
        --outdir . \\
        --max-runs ${maxRuns} \\
        ${longArg} ${emailArg} ${keyArg}
    test -s Runlist.tsv || { echo "Runlist.tsv is empty for '${species}'" >&2; exit 2; }
    """

    stub:
    """
    printf '@Run_acc\\ttotal_spots\\ttotal_bases\\tavg_len\\tbool:paired\\tcolor_space\\nSRR000001\\t100000\\t10000000\\t100.0\\t0\\t0\\n' > Runlist.tsv
    """
}


process VARUS_INDEX {
    tag { species }
    publishDir { "${params.outdir}/${species.replaceAll(' ', '_')}/varus" }, mode: 'copy', overwrite: true
    cpus { params.varus_index_cpus ?: 8 }
    // minimap2 (--longreads) builds parts of <= 8 Gbp; hisat2-build all at once
    memory { memGb(Math.max(8.0, 9 * (params.longreads ? Math.min(genomeGb(genome), 8.0) : genomeGb(genome)) + 4), task.attempt) }
    time { hours(Math.max(4.0, genomeGb(genome)), task.attempt) }
    errorStrategy { task.exitStatus in 137..140 && task.attempt <= 2 ? 'retry' : 'terminate' }

    input:
        tuple val(species), path(genome), path(runlist), val(extra)

    output:
        tuple val(species), path(genome), path(runlist),
              path("genome_index"), val(extra)

    script:
    def longArg = params.longreads ? '--longreads' : ''
    def prefix  = params.longreads ? 'mm2idx'      : 'hisatidx'
    """
    set -euo pipefail
    mkdir -p genome_index
    varus index ${genome} \\
        --outdir genome_index \\
        --threads ${task.cpus} \\
        --prefix ${prefix} \\
        ${longArg}
    """

    stub:
    """
    mkdir -p genome_index
    touch genome_index/hisatidx.1.ht2 genome_index/mm2idx.mmi
    """
}


process VARUS_LOGAN {
    // Pre-screen (on by default): align each candidate run's Logan contigs
    // (public S3, no credentials) to the genome, drop foreign runs, rank the
    // rest by tile coverage and seed the splice-site DB. --varus_logan false
    // skips it.
    tag { species }
    publishDir { "${params.outdir}/${species.replaceAll(' ', '_')}/varus" }, mode: 'copy', overwrite: true
    cpus { params.varus_logan_cpus ?: 8 }
    // enough to build the minimap2 index in one part (~8 bytes/base, the
    // fast mode; `varus logan` falls back to 8 Gbp parts in less memory)
    memory { memGb(Math.max(32.0, 9 * genomeGb(genome) + 16), task.attempt) }
    time { hours(Math.max(12.0, 2 * genomeGb(genome)), task.attempt) }
    errorStrategy { task.exitStatus in 137..140 && task.attempt <= 2 ? 'retry' : 'terminate' }

    input:
        tuple val(species), path(genome), path(runlist),
              path(index_dir), val(extra)

    output:
        tuple val(species), path(genome), path(runlist),
              path(index_dir), path("logan"), path("Runlist.logan.tsv"), val(extra)

    script:
    def maxCand   = params.varus_logan_max_candidates ?: 500
    def selectTop = params.varus_logan_select_top     ?: 50
    def mmiArg    = params.longreads ? "--mmi ${index_dir}/mm2idx.mmi" : ''
    """
    set -euo pipefail
    set +e
    varus logan ${genome} \\
        --runlist ${runlist} \\
        --outdir . \\
        --threads ${task.cpus} \\
        --max-candidates ${maxCand} \\
        --select-top ${selectTop} \\
        ${mmiArg}
    rc=\$?
    set -e
    # exit 3 (nothing accepted) and 4 (Logan unreachable) are not fatal.
    # 3: Runlist.logan.tsv still lists the runs Logan could not screen (too
    #    new for the last rebuild), which are often the only usable ones;
    #    keep them, but drop the ranking so VARUS_RUN applies no Logan prior.
    #    Fall back to the full runlist only if no run is left.
    # 4: nothing was written; VARUS_RUN gets the full runlist.
    case "\$rc" in
      0) ;;
      3) rm -f logan/LoganRanking.tsv
         if grep -q '^[^@]' Runlist.logan.tsv 2>/dev/null; then
           echo "varus logan: no run accepted for '${species}'; using the runs Logan could not screen" >&2
         else
           echo "varus logan: no run accepted for '${species}'; using the full runlist" >&2
           rm -f Runlist.logan.tsv
         fi ;;
      4) echo "varus logan: Logan S3 unreachable for '${species}'; using the full runlist" >&2 ;;
      *) exit "\$rc" ;;
    esac
    test -d logan || mkdir -p logan
    test -f Runlist.logan.tsv || cp ${runlist} Runlist.logan.tsv
    """

    stub:
    """
    mkdir -p logan
    cp ${runlist} Runlist.logan.tsv
    """
}


process VARUS_RUN {
    tag { species }
    publishDir { "${params.outdir}/${species.replaceAll(' ', '_')}/varus" }, mode: 'copy', overwrite: true
    cpus { params.varus_run_cpus ?: 16 }
    // aligner index (HISAT2 ~1.8 x the FASTA; minimap2 one part of <= 8 Gbp)
    // plus samtools sort buffers and the main process
    memory { memGb(Math.max(24.0, (params.longreads ? 4 * Math.min(genomeGb(genome), 8.0) : 2 * genomeGb(genome)) + 16), task.attempt) }

    input:
        tuple val(species), path(genome), path(runlist),
              path(index_dir), path(logan_dir), path(logan_runlist), val(extra)

    output:
        tuple val(species), path(genome), path("VARUS.bam"), val(extra), emit: bam
        path "introns.gff",       optional: true,                       emit: introns
        path "Coverage.csv",      optional: true,                       emit: coverage
        path "RunStatistics.csv", optional: true,                       emit: stats
        path "BatchTimings.tsv",  optional: true,                       emit: timings
        // archive these two with the genome: `varus replay` rebuilds VARUS.bam
        path "VARUS.manifest.tsv",                                      emit: manifest
        path "VARUS.splicedb.log.gz",                                   emit: splicedb_log
        path "runtime.varus.txt",                                       emit: runtime

    script:
    def maxBatches  = params.varus_max_batches ?: 1000
    // Long-read SRA runs have far fewer spots; default the batch-size lower
    // when the user hasn't overridden it.
    def defaultBatchSize = params.longreads ? 2000 : 50000
    def batchSize   = params.varus_batch_size  ?: defaultBatchSize
    def tileSize    = params.varus_tile_size   ?: 5000
    def minUniqPct  = params.varus_min_uniq_pct ?: 5.0
    def seed        = params.varus_seed         ?: 1
    def bootstrap   = params.varus_bootstrap_all ? '--bootstrap-all' : ''
    def profitCond  = params.varus_profit_condition ? '--profit-condition' : ''
    def parallelDl  = params.varus_parallel_downloads ?: 6
    def mergeEvery  = params.varus_merge_every != null ? params.varus_merge_every : 100
    def longArgs    = params.longreads ? '--longreads' : ''
    def indexPath   = params.longreads ? "${index_dir}/mm2idx.mmi" : "${index_dir}/hisatidx"
    def useLogan    = params.varus_logan ? true : false
    def loganTop    = params.varus_logan_top ?: 0
    """
    set -euo pipefail
    RUNLIST=${runlist}
    # The pre-screen ran (or was switched off) in VARUS_LOGAN, so `varus run`
    # must not start its own: --no-logan unless a ranking is handed over.
    # Runlist.logan.tsv is used whenever it has data rows; the prior
    # (--logan-dir) only when the pre-screen produced a ranking (exit 0).
    LOGAN_ARGS="--no-logan"
    if [ "${useLogan}" = "true" ] && grep -q '^[^@]' ${logan_runlist} 2>/dev/null; then
        RUNLIST=${logan_runlist}
        if [ -f ${logan_dir}/LoganRanking.tsv ]; then
            LOGAN_ARGS="--logan-dir ${logan_dir} --logan-top ${loganTop}"
        fi
    fi
    set +e
    /usr/bin/time -p -o runtime.varus.txt \\
      varus run '${species}' ${genome} \\
        --runlist \$RUNLIST \\
        --index ${indexPath} \\
        --outdir . \\
        --batch-size ${batchSize} \\
        --max-batches ${maxBatches} \\
        --tile-size ${tileSize} \\
        --min-uniq-pct ${minUniqPct} \\
        --threads ${task.cpus} \\
        --seed ${seed} \\
        --parallel-downloads ${parallelDl} \\
        --merge-every ${mergeEvery} \\
        ${bootstrap} ${profitCond} ${longArgs} \$LOGAN_ARGS
    rc=\$?
    set -e
    if [ "\$rc" = "3" ]; then
        echo "VARUS: no batch passed the quality gate for '${species}' (no usable RNA-seq)" >&2
        exit 3
    fi
    [ "\$rc" = "0" ] || exit "\$rc"
    test -s VARUS.bam || { echo "VARUS run produced no BAM" >&2; exit 2; }
    """

    stub:
    """
    touch VARUS.bam runtime.varus.txt introns.gff Coverage.csv RunStatistics.csv
    touch VARUS.manifest.tsv VARUS.splicedb.log.gz
    """
}

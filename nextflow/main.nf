#!/usr/bin/env nextflow
/*
 * Standalone VARUS v2 pipeline.
 *
 *   Input : a CSV with columns species,genome  (paths to per-species genome FASTAs)
 *   Output: per species  →  VARUS.bam, introns.gff, Coverage.csv, RunStatistics.csv
 *
 * Example:
 *   nextflow run nextflow/main.nf \
 *     --species_csv  examples/species.csv \
 *     --outdir       results \
 *     --ncbi_email   you@host
 */

nextflow.enable.dsl = 2

include { VARUS_RUNLIST; VARUS_INDEX; VARUS_LOGAN; VARUS_RUN } from './varus.nf'


// ---------------------------- params ----------------------------

params.species_csv          = params.species_csv          ?: null
params.outdir               = params.outdir               ?: 'results'

// VARUS run-time hyperparameters (all optional)
params.varus_max_batches    = (params.containsKey('varus_max_batches')   && params.varus_max_batches   != null ? params.varus_max_batches   : 1000) as int
params.varus_batch_size     = (params.containsKey('varus_batch_size')    && params.varus_batch_size    != null ? params.varus_batch_size    : 50000) as int
params.varus_tile_size      = (params.containsKey('varus_tile_size')     && params.varus_tile_size     != null ? params.varus_tile_size     : 5000) as int
params.varus_min_uniq_pct   = (params.containsKey('varus_min_uniq_pct')  && params.varus_min_uniq_pct  != null ? params.varus_min_uniq_pct  : 5.0) as double
params.varus_max_runs       = (params.containsKey('varus_max_runs')      && params.varus_max_runs      != null ? params.varus_max_runs      : 0) as int
params.varus_seed           = (params.containsKey('varus_seed')          && params.varus_seed          != null ? params.varus_seed          : 1) as int
params.varus_bootstrap_all  = (params.containsKey('varus_bootstrap_all') ? params.varus_bootstrap_all : false) as boolean
params.varus_profit_condition = (params.containsKey('varus_profit_condition') ? params.varus_profit_condition : false) as boolean
params.varus_index_cpus     = (params.containsKey('varus_index_cpus')    && params.varus_index_cpus    != null ? params.varus_index_cpus    : 8) as int
params.varus_run_cpus       = (params.containsKey('varus_run_cpus')      && params.varus_run_cpus      != null ? params.varus_run_cpus      : 16) as int

// Speed knobs (v2): concurrent batch downloads, .sra prefetch, rolling merge.
params.varus_parallel_downloads = (params.containsKey('varus_parallel_downloads') && params.varus_parallel_downloads != null ? params.varus_parallel_downloads : 6) as int
params.varus_prefetch       = (params.containsKey('varus_prefetch') ? params.varus_prefetch : false) as boolean
params.varus_merge_every    = (params.containsKey('varus_merge_every') && params.varus_merge_every != null ? params.varus_merge_every : 100) as int

// Logan pre-screen (v2): on by default (--varus_logan false to skip); needs minimap2.
params.varus_logan          = (params.containsKey('varus_logan') ? params.varus_logan : true) as boolean
params.varus_logan_cpus     = (params.containsKey('varus_logan_cpus') && params.varus_logan_cpus != null ? params.varus_logan_cpus : 8) as int
params.varus_logan_max_candidates = (params.containsKey('varus_logan_max_candidates') && params.varus_logan_max_candidates != null ? params.varus_logan_max_candidates : 500) as int
params.varus_logan_select_top = (params.containsKey('varus_logan_select_top') && params.varus_logan_select_top != null ? params.varus_logan_select_top : 50) as int
params.varus_logan_top      = (params.containsKey('varus_logan_top') && params.varus_logan_top != null ? params.varus_logan_top : 0) as int

// Long-read mode: align with minimap2, restrict the SRA query to PacBio/ONT.
// Implies a different splice-DB format and a smaller default --batch-size.
// The minimap2 preset is auto-selected per run from SRA platform metadata.
params.longreads            = (params.containsKey('longreads') ? params.longreads : false) as boolean

params.ncbi_email           = params.ncbi_email   ?: null
params.ncbi_api_key         = params.ncbi_api_key ?: null


def die(msg) { log.error msg; System.exit(1) }

if (!params.species_csv) die("Missing --species_csv")


// ---------------------------- workflow ----------------------------

workflow {

    ch_input = Channel.fromPath(params.species_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.species, file(row.genome), [:]) }
        // (species, genome_path, extra) — `extra` is a placeholder so callers
        // wiring this into a larger workflow can carry context downstream.

    runlist_out = VARUS_RUNLIST(ch_input)
    index_out   = VARUS_INDEX(runlist_out)
    if (params.varus_logan) {
        logan_out = VARUS_LOGAN(index_out)
    } else {
        // Same tuple shape without a pre-screen: empty logan dir, same runlist.
        logan_out = index_out.map { species, genome, runlist, index_dir, extra ->
            tuple(species, genome, runlist, index_dir, file("$projectDir/NO_LOGAN", type: 'dir'), runlist, extra)
        }
    }
    bam_out     = VARUS_RUN(logan_out).bam
}

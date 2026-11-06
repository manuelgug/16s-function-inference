#!/usr/bin/env nextflow
// GTDB genome pipeline: annotate genomes (truth for the real-world test, DESIGN §B.2) and Aim 2 (paused).
// Start it with ./run.sh test|pilot|full, which sets every parameter. The benchmark is benchmark.nf.
// Stages: select genomes -> per-batch fetch / barrnap / genes+16S / eggNOG -> targets
//         and SINA alignment -> allele features -> QC report -> Aim 2.

include { SELECT_GENOMES        } from './modules/select_genomes'
include { FETCH_GENOMES         } from './modules/fetch_genomes'
include { BARRNAP               } from './modules/barrnap'
include { PROCESS_GENOMES       } from './modules/process_genomes'
include { EGGNOG                } from './modules/eggnog'
include { BUILD_TARGETS         } from './modules/build_targets'
include { SINA_INDEX; ALIGN_16S } from './modules/sina'
include { BUILD_FEATURES        } from './modules/build_features'
include { QC_REPORT             } from './modules/qc_report'
include { AIM2_PREP; AIM2_PYSEER; AIM2_REPORT } from './modules/aim2'

workflow {
    def required = ['gtdb_release', 'silva_release', 'min_completeness', 'max_contamination', 'min_ssu_len',
                    'full_n', 'seed', 'min_per_phylum', 'min_col_coverage', 'min_allele_freq', 'primers',
                    'pcr_error_rate', 'ko_min_prevalence', 'ko_max_prevalence',
                    'batch_size', 'max_cpus', 'max_memory', 'db_dir', 'outdir']
    def missing = required.findAll { params[it] == null }
    if (missing) {
        error "Missing parameters: ${missing.join(', ')}. Start the pipeline with ./run.sh, which sets them."
    }
    def db   = params.db_dir
    def gtdb = "${db}/gtdb_r${params.gtdb_release}"

    // --- genome list: committed list (test, pilot) or filtered from GTDB (full) ---
    genome_list = params.genome_list
        ? channel.fromPath(params.genome_list, checkIfExists: true)
        : SELECT_GENOMES(file("${gtdb}/bac120_metadata_r${params.gtdb_release}.tsv.gz", checkIfExists: true),
                         file("${gtdb}/ar53_metadata_r${params.gtdb_release}.tsv.gz", checkIfExists: true))

    batches = genome_list
        .splitText(by: params.batch_size as int, keepHeader: true, file: true)
        .map { f -> tuple(f.baseName, f) }

    // --- per batch ---
    FETCH_GENOMES(batches)
    BARRNAP(FETCH_GENOMES.out.genomes)
    PROCESS_GENOMES(FETCH_GENOMES.out.genomes.join(BARRNAP.out))
    EGGNOG(PROCESS_GENOMES.out.faa.filter { _batch, faa -> faa.size() > 0 },
           file("${db}/eggnog_5.0.2", type: 'dir', checkIfExists: true))

    // --- whole set ---
    BUILD_TARGETS(EGGNOG.out.collect())

    ssu   = PROCESS_GENOMES.out.ssu.collectFile(name: 'ssu.fna', sort: { f -> f.name })
    ecoli = file("${projectDir}/config/ecoli_rrsB_16S.fna", checkIfExists: true)
    SINA_INDEX(file("${db}/silva_${params.silva_release}/SILVA_${params.silva_release}_SSURef_NR99_opt.arb", checkIfExists: true),
               ecoli)
    ALIGN_16S(ssu, ecoli, SINA_INDEX.out)
    BUILD_FEATURES(ALIGN_16S.out.aln, ssu, file(params.primers, checkIfExists: true))

    QC_REPORT(genome_list,
              FETCH_GENOMES.out.status.collectFile(name: 'fetch_status.tsv', keepHeader: true, sort: { f -> f.name }),
              PROCESS_GENOMES.out.stats.collectFile(name: 'genome_stats.tsv', keepHeader: true, sort: { f -> f.name }),
              BUILD_TARGETS.out, BUILD_FEATURES.out, ALIGN_16S.out.meta)

    // --- Aim 2: association beyond phylogeny (pyseer LMM), bacteria only ---
    AIM2_PREP(file("${gtdb}/bac120_r${params.gtdb_release}.tree", checkIfExists: true), genome_list,
              QC_REPORT.out, BUILD_FEATURES.out, BUILD_TARGETS.out)
    AIM2_PYSEER(AIM2_PREP.out.chunks.flatten(), AIM2_PREP.out.dir.first())
    AIM2_REPORT(AIM2_PYSEER.out.collect(), AIM2_PREP.out.dir, BUILD_FEATURES.out)
}

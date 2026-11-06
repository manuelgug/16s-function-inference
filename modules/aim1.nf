// Cross-validation folds at each holdout level, over usable genomes only.
process MAKE_SPLITS {
    publishDir params.outdir, mode: 'copy'

    input:
    path genomes
    path qc

    output:
    path 'splits.tsv'

    script:
    """
    make_splits.py --genomes ${genomes} --usable ${qc}/genome_qc.tsv --levels ${params.holdout_levels} \\
        --folds ${params.cv_folds} --seed ${params.seed} --out splits.tsv
    """
}

// Aim 1 models for one input arm (full or an amplicon) and one holdout level.
process AIM1_PREDICT {
    tag "${arm}/${level}"
    publishDir "${params.outdir}/aim1", mode: 'copy'

    input:
    tuple val(arm), val(level)
    path features
    path targets
    path splits

    output:
    path "${arm}_${level}"

    script:
    """
    aim1_predict.py --features ${features} --targets ${targets} --splits ${splits} --arm ${arm} --level ${level} \\
        --ko-min-prev ${params.ko_min_prevalence} --ko-max-prev ${params.ko_max_prevalence} \\
        --threads ${task.cpus} --seed ${params.seed} --outdir ${arm}_${level}
    """
}

// Tables, figures and markdown report over all Aim 1 tasks; also copied into the repo's reports/ folder.
process AIM1_REPORT {
    publishDir params.outdir, mode: 'copy'
    publishDir "${projectDir}/reports/${workflow.profile.tokenize(',')[0]}", mode: 'copy'

    input:
    path runs
    path patristic
    path splits
    path picrust2

    output:
    path 'aim1_report'

    script:
    """
    aim1_report.py --runs ${runs} --outdir aim1_report --seed ${params.seed} \
        --patristic ${patristic} --splits ${splits} --picrust2 ${picrust2}
    render_html.py aim1_report/aim1_report.md
    """
}

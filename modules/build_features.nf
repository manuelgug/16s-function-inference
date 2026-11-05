// Allele-at-position features for full-length 16S and each in-silico amplicon.
process BUILD_FEATURES {
    publishDir params.outdir, mode: 'copy'

    input:
    path aln
    path ssu
    path primers

    output:
    path 'features'

    script:
    """
    build_features.py --aln ${aln} --ssu ${ssu} --primers ${primers} --ref-id ecoli_rrsB \\
        --min-col-coverage ${params.min_col_coverage} --min-allele-freq ${params.min_allele_freq} \\
        --pcr-error-rate ${params.pcr_error_rate} --outdir features
    """
}

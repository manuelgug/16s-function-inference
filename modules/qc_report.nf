// Per-genome yield table and a short markdown summary; also copied into the repo's reports/ folder.
process QC_REPORT {
    publishDir params.outdir, mode: 'copy'
    publishDir "${projectDir}/reports/${workflow.profile.tokenize(',')[0]}", mode: 'copy'

    input:
    path genomes
    path fetch_status
    path stats
    path targets
    path features
    path sina_meta

    output:
    path 'qc'

    script:
    """
    qc_report.py --genomes ${genomes} --fetch ${fetch_status} --stats ${stats} \\
        --targets ${targets} --features ${features} --sina ${sina_meta} --outdir qc
    render_html.py qc/qc_report.md
    """
}

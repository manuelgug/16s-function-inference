// Genome x function matrices: KO presence, KEGG module completeness, COG category profile.
process BUILD_TARGETS {
    publishDir params.outdir, mode: 'copy'

    input:
    path annotations

    output:
    path 'targets'

    script:
    """
    build_targets.py --annotations ${annotations} --outdir targets --cpus ${task.cpus}
    """
}

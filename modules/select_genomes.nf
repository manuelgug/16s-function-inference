// Full scale only: filter GTDB representatives and optionally draw a stratified sample (FULL_N).
process SELECT_GENOMES {
    publishDir params.outdir, mode: 'copy'

    input:
    path bac_meta
    path arc_meta

    output:
    path 'genomes.tsv'

    script:
    """
    select_genomes.py --bac ${bac_meta} --arc ${arc_meta} \\
        --min-completeness ${params.min_completeness} --max-contamination ${params.max_contamination} \\
        --min-ssu-len ${params.min_ssu_len} --n ${params.full_n} \\
        --min-per-phylum ${params.min_per_phylum} --seed ${params.seed} --out genomes.tsv
    """
}

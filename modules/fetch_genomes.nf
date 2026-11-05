// Download one batch of genome assemblies from NCBI. Failed genomes are reported, not fatal.
process FETCH_GENOMES {
    tag "${batch}"

    input:
    tuple val(batch), path(genomes)

    output:
    tuple val(batch), path(genomes), path('fna'), emit: genomes
    path "${batch}.fetch_status.tsv",              emit: status

    script:
    """
    fetch_genomes.py ${genomes} fna --threads 4
    mv fetch_status.tsv ${batch}.fetch_status.tsv
    """
}

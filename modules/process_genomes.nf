// Gene calling (pyrodigal) and representative 16S selection for one batch.
process PROCESS_GENOMES {
    tag "${batch}"

    input:
    tuple val(batch), path(genomes), path(fna), path(gff)

    output:
    tuple val(batch), path("${batch}.faa"), emit: faa
    path "${batch}.ssu.fna",                emit: ssu
    path "${batch}.ssu_copies.tsv",         emit: copies
    path "${batch}.stats.tsv",              emit: stats

    script:
    """
    process_genomes.py --genomes ${genomes} --fna-dir ${fna} --gff-dir ${gff} --prefix ${batch} \\
        --min-ssu-len ${params.min_ssu_len} --threads ${task.cpus}
    """
}

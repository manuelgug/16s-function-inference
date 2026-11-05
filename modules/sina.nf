// One-off: copy the SILVA ARB file and let SINA build its k-mer index (.sidx) next to it.
// storeDir keeps the result in the database folder, so this runs once per machine.
process SINA_INDEX {
    storeDir "${params.db_dir}/silva_${params.silva_release}/sina_index"

    input:
    path arb
    path probe

    output:
    tuple path('silva.arb'), path('silva.sidx')

    script:
    """
    cp -L ${arb} silva.arb
    sina -i ${probe} -o probe.aligned.fasta --db silva.arb -p ${task.cpus}
    """
}

// Align every representative 16S plus the E. coli reference to the SILVA reference alignment.
process ALIGN_16S {
    publishDir params.outdir, mode: 'copy', pattern: 'alignment/*'

    input:
    path ssu
    path ecoli
    tuple path(arb), path(sidx)

    output:
    path 'alignment/aligned.fasta', emit: aln
    path 'alignment/aligned.csv',   emit: meta

    script:
    """
    mkdir alignment
    cat ${ssu} ${ecoli} > input.fna
    sina -i input.fna -o alignment/aligned.fasta --db ${arb} -p ${task.cpus} \\
        --meta-fmt csv --fields align_ident_slv,align_quality_slv,align_cutoff_head_slv,align_cutoff_tail_slv
    """
}

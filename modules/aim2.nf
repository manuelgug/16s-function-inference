// Aim 2 inputs: pruned GTDB tree, kinship, allele genotypes, phenotype chunks (bacteria only, D019).
process AIM2_PREP {
    input:
    path tree
    path genomes
    path qc
    path features
    path targets

    output:
    path 'aim2_input', emit: dir
    path 'aim2_input/pheno_*.tsv', emit: chunks

    script:
    """
    aim2_prep.py --tree ${tree} --genomes ${genomes} --qc ${qc}/genome_qc.tsv --features ${features} \\
        --targets ${targets} --ko-min-prev ${params.ko_min_prevalence} --ko-max-prev ${params.ko_max_prevalence} \\
        --chunks 32 --outdir aim2_input
    """
}

// pyseer LMM, every phenotype in one chunk file against every 16S allele.
process AIM2_PYSEER {
    tag "${pheno.baseName}"

    input:
    path pheno
    path input_dir

    output:
    path "${pheno.baseName}.pyseer.tsv.gz"

    script:
    def cont = pheno.name.contains('continuous') ? '--continuous' : ''
    """
    for ph in \$(head -1 ${pheno} | cut -f2-); do
        pyseer --lmm ${cont} --phenotypes ${pheno} --phenotype-column "\$ph" --pres ${input_dir}/genotypes.Rtab \\
            --similarity ${input_dir}/kinship.tsv --cpu 1 2>> pyseer.log \\
          | awk -v p="\$ph" 'BEGIN{OFS="\\t"} NR>1 {print p, \$1, \$2, \$4, \$5, \$8}'
    done | gzip > ${pheno.baseName}.pyseer.tsv.gz
    """
}

process AIM2_REPORT {
    publishDir params.outdir, mode: 'copy'
    publishDir "${projectDir}/reports/${workflow.profile.tokenize(',')[0]}", mode: 'copy'

    input:
    path results
    path input_dir
    path features

    output:
    path 'aim2_report'

    script:
    """
    aim2_report.py --results ${results} --input ${input_dir} --features ${features} --outdir aim2_report
    """
}

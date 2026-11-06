// M2: PICRUSt2 out of the box on every representative 16S (external reference only, DESIGN §4.1).
process PICRUST2 {
    publishDir params.outdir, mode: 'copy'   // results/<profile>/picrust2_out/: raw PICRUSt2 predictions + NSTI
    input:
    path ssu

    output:
    tuple path('picrust2_out/combined_KO_predicted.tsv.gz'), path('picrust2_out/combined_marker_predicted_and_nsti.tsv.gz'),
          path('picrust2_out/reference_genomes.txt'), emit: pred
    path 'picrust2_out/ko_info.tsv.gz', emit: ko_names

    script:
    """
    export HOME=\$PWD   # Apptainer --no-home leaves HOME read-only
    ( printf 'sequence\\tsample1\\n'; grep '>' ${ssu} | sed 's/^>//; s/\$/\\t1/' ) > table.tsv
    picrust2_pipeline.py -s ${ssu} -i table.tsv -o picrust2_out -p ${task.cpus} --skip_norm --in_traits KO --no_pathways
    # genomes in PICRUSt2's own reference (bacteria + archaea trees): a test genome listed here is leakage
    ref=\$(python -c 'import picrust2, os; print(os.path.dirname(picrust2.__file__))')/default_files
    grep -h '>' \$ref/bacteria/bac_ref/bac_ref.fna \$ref/archaea/arc_ref/arc_ref.fna | sed 's/^>//; s/ .*//' > picrust2_out/reference_genomes.txt
    cp \$ref/description_mapfiles/ko_info.tsv.gz picrust2_out/   # KO names, for readable tables
    """
}

process PICRUST2_EVAL {
    input:
    tuple path(pred), path(nsti), path(refs)
    path targets

    output:
    path 'picrust2_metrics.tsv'

    script:
    """
    picrust2_eval.py --pred ${pred} --nsti ${nsti} --refs ${refs} --targets ${targets} \\
        --ko-min-prev ${params.ko_min_prevalence} --ko-max-prev ${params.ko_max_prevalence} --out picrust2_metrics.tsv
    """
}

// Head-to-head with PICRUSt2 on the same genomes (answers: who wins, where, by how much; ordination).
process PICRUST2_COMPARE {
    publishDir params.outdir, mode: 'copy'
    publishDir "${projectDir}/reports/${workflow.profile.tokenize(',')[0]}", mode: 'copy',
        saveAs: { f -> f.endsWith('.gz') ? null : f }   // the large call tables stay in results/

    input:
    tuple path(pred), path(nsti), path(refs)
    path ko_names
    path targets
    path aim1_full_random
    path genomes

    output:
    path 'picrust2_comparison/*'   // glob, so saveAs sees each file name

    script:
    """
    picrust2_compare.py --pred ${pred} --nsti ${nsti} --refs ${refs} --ko-names ${ko_names} --targets ${targets} \\
        --aim1 ${aim1_full_random} --genomes ${genomes} --cpus ${task.cpus} --seed ${params.seed} \\
        --outdir picrust2_comparison
    render_html.py picrust2_comparison/picrust2_comparison.md
    """
}

// GTDB-tree (patristic) distances between usable genomes, for distance-to-training curves.
process PATRISTIC {
    input:
    path bac_tree
    path arc_tree
    path qc

    output:
    path 'patristic.parquet'

    script:
    """
    patristic.py --bac-tree ${bac_tree} --arc-tree ${arc_tree} --qc ${qc}/genome_qc.tsv --out patristic.parquet
    """
}

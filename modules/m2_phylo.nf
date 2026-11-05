// M2: PICRUSt2 out of the box on every representative 16S (external reference only, DESIGN §4.1).
process PICRUST2 {
    input:
    path ssu

    output:
    tuple path('picrust2_out/combined_KO_predicted.tsv.gz'), path('picrust2_out/combined_marker_predicted_and_nsti.tsv.gz')

    script:
    """
    export HOME=\$PWD   # Apptainer --no-home leaves HOME read-only
    ( printf 'sequence\\tsample1\\n'; grep '>' ${ssu} | sed 's/^>//; s/\$/\\t1/' ) > table.tsv
    picrust2_pipeline.py -s ${ssu} -i table.tsv -o picrust2_out -p ${task.cpus} --skip_norm --in_traits KO --no_pathways
    """
}

process PICRUST2_EVAL {
    input:
    tuple path(pred), path(nsti)
    path targets

    output:
    path 'picrust2_metrics.tsv'

    script:
    """
    picrust2_eval.py --pred ${pred} --nsti ${nsti} --targets ${targets} \\
        --ko-min-prev ${params.ko_min_prevalence} --ko-max-prev ${params.ko_max_prevalence} --out picrust2_metrics.tsv
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

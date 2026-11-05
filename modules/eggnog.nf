// Functional annotation of one batch's proteins with eggNOG-mapper (DIAMOND mode, default sensitivity).
process EGGNOG {
    tag "${batch}"

    input:
    tuple val(batch), path(faa)
    path eggnog_db

    output:
    path "${batch}.emapper.annotations"

    script:
    """
    emapper.py -i ${faa} --itype proteins -m diamond --data_dir ${eggnog_db} --cpu ${task.cpus} \\
        --output ${batch} --output_dir . --temp_dir . --override
    """
}

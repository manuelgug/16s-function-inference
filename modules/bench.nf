// Benchmark (DESIGN §B, BENCHMARK.md): s16fun vs PICRUSt2 and baselines on PICRUSt2's own reference.

// Nextflow's cache key ignores the content of bin/ scripts: each process embeds a hash of the scripts it runs,
// so changing code re-runs exactly the tasks that use it (no stale results, D091).
def code(List names) { names.collect { n -> file("${projectDir}/bin/${n}").text }.join('').md5() }

// PICRUSt2's shipped reference (PICRUSt2-SC), copied out of its container.
process PSC_EXPORT {
    output:
    path 'psc'

    script:
    """
    ref=\$(python -c 'import picrust2, os; print(os.path.dirname(picrust2.__file__))')/default_files
    mkdir -p psc/bacteria psc/archaea
    for d in bacteria archaea; do
        cp -r \$ref/\$d/\${d:0:3}_ref \$ref/\$d/ko.txt.gz \$ref/\$d/16S.txt.gz \$ref/\$d/\${d}_metadata.csv.gz psc/\$d/
    done
    cp \$ref/description_mapfiles/ko_info.tsv.gz psc/
    """
}

process BENCH_EXTRACT {
    publishDir params.outdir, mode: 'copy', saveAs: { f -> f == 'bench_data/genomes.tsv' ? f : null }
    input:
    path psc

    output:
    path 'bench_data'

    script:
    """
    # code ${code(['bench_extract.py'])}
    bench_extract.py --psc ${psc} --outdir bench_data --min-ko-genomes ${params.bench_min_ko_genomes} \\
        --subsample ${params.bench_subsample} --seed ${params.seed}
    """
}

process BENCH_SPLITS {
    publishDir params.outdir, mode: 'copy'
    input:
    path data

    output:
    path 'splits.tsv'

    script:
    """
    # code ${code(['make_splits.py'])}
    make_splits.py --genomes ${data}/genomes.tsv --levels ${params.holdout_levels} --folds ${params.cv_folds} \\
        --seed ${params.seed} --dev-fraction ${params.bench_dev_fraction} --out splits.tsv
    """
}

// Training set + query sequences of one level x fold x arm (identical for every method).
process BENCH_PREPARE {
    tag "${level}/${fold}/${arm}"
    input:
    tuple val(level), val(fold), val(arm)
    path data
    path splits
    path primers

    output:
    tuple val(level), val(fold), val(arm), path('task')

    script:
    """
    # code ${code(['bench_prepare.py', 'benchlib.py', 'build_features.py'])}
    bench_prepare.py --data ${data} --splits ${splits} --level ${level} --fold ${fold} --arm ${arm} \\
        --primers ${primers} --pcr-error-rate ${params.pcr_error_rate} --phase ${params.bench_phase} --outdir task
    """
}

// PICRUSt2 reference files restricted to one fold's training genomes (shared by all arms).
process BENCH_P2_REF {
    tag "${level}/${fold}"
    input:
    tuple val(level), val(fold), path(tdir)
    path psc

    output:
    tuple val(level), val(fold), path('ref')

    script:
    """
    # code ${code(['bench_picrust2_ref.py', 'bench_extract.py'])}
    bench_picrust2_ref.py --psc ${psc} --train ${tdir}/train.txt --outdir ref
    """
}

// PICRUSt2 2.6.4 on one task: placement once, then every HSP setting on the placed trees.
process PICRUST2_FOLD {
    tag "${level}/${fold}/${arm}"
    input:
    tuple val(level), val(fold), val(arm), path(tdir), path(ref)

    output:
    tuple val(level), val(fold), val(arm), path('p2')

    script:
    def methods = params.bench_p2_methods.tokenize(',')
    def others = methods.findAll { m -> m != 'mp' }.join(' ')
    """
    export HOME=\$PWD   # Apptainer --no-home leaves HOME read-only
    cp -rL ${ref} r && mkdir p2
    for d in bac arc; do
        hmmbuild --dna --informat afa --cpu ${task.cpus} r/\${d}_ref/\${d}_ref.hmm r/\${d}_ref/\${d}_ref.fna > /dev/null
    done
    if [ ! -s ${tdir}/query.fna ]; then exit 0; fi
    ( printf 'sequence\\tsample1\\n'; grep '>' ${tdir}/query.fna | sed 's/^>//; s/\$/\\t1/' ) > table.tsv
    picrust2_pipeline.py -s ${tdir}/query.fna -i table.tsv -o out -p ${task.cpus} --skip_norm --no_pathways \\
        -r1 r/bac_ref -r2 r/arc_ref --custom_trait_tables_ref1 r/bac_ko.tsv.gz --custom_trait_tables_ref2 r/arc_ko.tsv.gz \\
        --marker_gene_table_ref1 r/bac_16S.tsv.gz --marker_gene_table_ref2 r/arc_16S.tsv.gz
    cp out/combined_marker_predicted_and_nsti.tsv.gz r/bac_constant.tsv r/arc_constant.tsv p2/
    mv p2/combined_marker_predicted_and_nsti.tsv.gz p2/marker_nsti.tsv.gz
    for i in 1 2; do
        d=\$([ \$i = 1 ] && echo bac || echo arc)
        [ -s out/ref\${i}_\${d}_ko.tsv_predicted.tsv.gz ] && cp out/ref\${i}_\${d}_ko.tsv_predicted.tsv.gz p2/mp_ref\${i}.tsv.gz
        for m in ${others}; do
            [ -s out/ref\${i}_reduced.tre ] || continue
            hsp.py --tree out/ref\${i}_reduced.tre --observed_trait_table r/\${d}_ko.tsv.gz --hsp_method \$m \\
                --edge_exponent 0.5 --seed 100 --processes ${task.cpus} --output p2/\${m}_ref\${i}.tsv.gz
        done
    done
    """
}

process BENCH_P2_PRED {
    tag "${level}/${fold}/${arm}"
    input:
    tuple val(level), val(fold), val(arm), path(tdir), path(p2)
    path data

    output:
    tuple val(level), val(fold), val(arm), path('pred_p2')

    script:
    """
    # code ${code(['bench_picrust2_pred.py', 'benchlib.py'])}
    bench_picrust2_pred.py --data ${data} --task ${tdir} --p2 ${p2} --methods ${params.bench_p2_methods} --outdir pred_p2
    """
}

// s16fun evidence + nearest-neighbour baselines (pass 1; dev also samples stacker training rows).
process BENCH_PREDICT {
    tag "${level}/${fold}/${arm}"
    input:
    tuple val(level), val(fold), val(arm), path(tdir)
    path data

    output:
    tuple val(level), val(fold), val(arm), path('pred1'), emit: pred
    path 'pred1/stacker_rows.npz', optional: true, emit: rows
    path 'pred1/ssu_rows.npz', optional: true, emit: ssu_rows

    script:
    """
    # code ${code(['bench_predict.py', 's16fun.py', 'benchlib.py'])}
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus}
    bench_predict.py --data ${data} --task ${tdir} --threads ${task.cpus} --seed ${params.seed} \\
        --sample-genomes ${params.bench_stacker_genomes} --outdir pred1
    """
}

process BENCH_MLP {
    tag "${level}/${fold}/${arm}"
    input:
    tuple val(level), val(fold), val(arm), path(tdir)
    path data

    output:
    tuple val(level), val(fold), val(arm), path('pred_mlp')

    script:
    """
    # code ${code(['bench_mlp.py', 'benchlib.py'])}
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus}
    bench_mlp.py --data ${data} --task ${tdir} --seed ${params.seed} --outdir pred_mlp
    """
}

// Stacker, cross-fitted by family on dev rows (dev phase only; confirm uses the frozen models).
process BENCH_STACKER {
    publishDir "${params.outdir}/frozen", mode: 'copy'
    input:
    path rows, stageAs: 'rows*/stacker_rows.npz'
    path ssu_rows, stageAs: 'ssu*/ssu_rows.npz'
    path data
    path splits

    output:
    path 'stackers'

    script:
    """
    # code ${code(['bench_stacker.py', 's16fun.py'])}
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus}
    bench_stacker.py --rows ${rows} --ssu-rows ${ssu_rows} --genomes ${data}/genomes.tsv --splits ${splits} --threads ${task.cpus} \\
        --seed ${params.seed} --outdir stackers
    """
}

// s16fun as shipped (pass 2) + ablations.
process BENCH_PREDICT_STACKED {
    tag "${level}/${fold}/${arm}"
    input:
    tuple val(level), val(fold), val(arm), path(tdir)
    path data
    path stackers

    output:
    tuple val(level), val(fold), val(arm), path('pred2')

    script:
    """
    # code ${code(['bench_predict.py', 's16fun.py', 'benchlib.py'])}
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus}
    bench_predict.py --data ${data} --task ${tdir} --stackers ${stackers} --threads ${task.cpus} --outdir pred2
    """
}

// One directory per task with every method's predictions, for the evaluation.
process BENCH_MERGE {
    tag "${level}/${fold}/${arm}"
    input:
    tuple val(level), val(fold), val(arm), path(tdir), path(p2), path(pred1), path(mlp), path(pred2)

    output:
    path "${level}_${fold}_${arm}"

    script:
    """
    mkdir ${level}_${fold}_${arm}
    cp -L ${tdir}/task.tsv ${tdir}/train.txt ${level}_${fold}_${arm}/
    cp -L ${p2}/* ${pred1}/pred_*.npz ${mlp}/pred_*.npz ${pred2}/pred_*.npz ${pred2}/novelty.tsv ${level}_${fold}_${arm}/
    """
}

// Simulated communities from genus-level held-out genomes (exact truth).
process BENCH_COMMUNITY {
    input:
    path tasks
    path data

    output:
    path 'community'

    script:
    """
    # code ${code(['bench_community.py', 'benchlib.py', 'bench_figures.py'])}
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus}
    bench_community.py --data ${data} --tasks ${tasks} --level genus --methods ${params.bench_p2_methods} \\
        --seed ${params.seed} --outdir community
    """
}

process BENCH_EVAL {
    publishDir "${projectDir}/reports/${workflow.profile.tokenize(',')[0]}", mode: 'copy',
        saveAs: { f -> f.endsWith('.parquet') || f.endsWith('.tsv.gz') ? null : f.replaceFirst(/^report/, params.bench_phase) }
    publishDir params.outdir, mode: 'copy'
    input:
    path tasks
    path data
    path splits
    path frozen
    path community
    path traces

    output:
    path 'report/*'

    script:
    def thr = params.bench_phase == 'confirm' ? "--thresholds ${frozen}/picrust2_thresholds.tsv" : ''
    """
    # code ${code(['bench_eval.py', 'bench_report.py', 'bench_figures.py', 'bench_ordination.py', 'benchlib.py', 'build_targets.py', 'render_html.py'])}
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus}
    bench_eval.py --data ${data} --splits ${splits} --tasks ${tasks} --phase ${params.bench_phase} ${thr} \\
        --stacker-importance ${frozen}/stacker_importance.tsv --community ${community} --trace ${traces} --boot ${params.bench_boot} \\
        --seed ${params.seed} --cpus ${task.cpus} --outdir report
    render_html.py report/benchmark_report.md
    """
}

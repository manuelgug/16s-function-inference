// Benchmark workflow (DESIGN §B, BENCHMARK.md). Started by ./run.sh bench_test | bench.
// Every method is scored on the same held-out genomes of PICRUSt2's own reference; see modules/bench.nf.

include { PSC_EXPORT; BENCH_EXTRACT; BENCH_SPLITS; BENCH_PREPARE; BENCH_P2_REF; PICRUST2_FOLD; BENCH_P2_PRED;
          BENCH_PREDICT; BENCH_MLP; BENCH_STACKER; BENCH_PREDICT_STACKED; BENCH_MERGE; BENCH_COMMUNITY; BENCH_EVAL } from './modules/bench'

workflow {
    def required = ['seed', 'holdout_levels', 'cv_folds', 'primers', 'pcr_error_rate', 'bench_phase', 'bench_arms',
                    'bench_p2_methods', 'bench_min_ko_genomes', 'bench_dev_fraction', 'bench_subsample',
                    'bench_stacker_genomes', 'bench_boot', 'bench_frozen', 'max_cpus', 'max_memory', 'outdir']
    def missing = required.findAll { params[it] == null }
    if (missing) {
        error "Missing parameters: ${missing.join(', ')}. Start the benchmark with ./run.sh, which sets them."
    }
    if (!(params.bench_phase in ['dev', 'confirm'])) {
        error "bench_phase must be dev or confirm"
    }

    PSC_EXPORT()
    data = BENCH_EXTRACT(PSC_EXPORT.out).first()
    splits = BENCH_SPLITS(data).first()
    primers = file(params.primers, checkIfExists: true)

    tasks = channel.fromList(params.holdout_levels.tokenize(','))
        .combine(channel.of(0..<(params.cv_folds as int)).flatten())
        .combine(channel.fromList(params.bench_arms.tokenize(',')))
    BENCH_PREPARE(tasks, data, splits, primers)
    prepared = BENCH_PREPARE.out

    // PICRUSt2: one reference per level x fold (the full-arm task carries the training list)
    refs = BENCH_P2_REF(prepared.filter { _l, _f, arm, _t -> arm == 'full' }.map { l, f, _a, t -> tuple(l, f, t) },
                        PSC_EXPORT.out.first())
    p2 = PICRUST2_FOLD(prepared.combine(refs, by: [0, 1]))
    p2pred = BENCH_P2_PRED(prepared.join(p2, by: [0, 1, 2]), data)

    // our method, baselines, k-mer net
    BENCH_PREDICT(prepared, data)
    mlp = BENCH_MLP(prepared, data)
    stackers = params.bench_phase == 'dev'
        ? BENCH_STACKER(BENCH_PREDICT.out.rows.collect(), BENCH_PREDICT.out.ssu_rows.collect(), data, splits).first()
        : channel.value(file(params.bench_frozen, type: 'dir', checkIfExists: true))
    pred2 = BENCH_PREDICT_STACKED(prepared, data, stackers)

    merged = BENCH_MERGE(prepared.join(p2pred, by: [0, 1, 2]).join(BENCH_PREDICT.out.pred, by: [0, 1, 2])
                         .join(mlp, by: [0, 1, 2]).join(pred2, by: [0, 1, 2]))
    traces = channel.fromPath("${params.outdir}/pipeline_info/trace_*.tsv").collect()
    community = BENCH_COMMUNITY(merged.collect(), data)
    BENCH_EVAL(merged.collect(), data, splits, stackers, community, traces)
}

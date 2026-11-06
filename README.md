# 16s-function-inference

Predict the functional repertoire (KEGG KOs and modules) of prokaryotes from their 16S rRNA gene, and benchmark
that prediction against PICRUSt2 and other 16S function predictors on identical data, with no leakage.

- **Design and plan:** [DESIGN.md](DESIGN.md) (benchmark: §B; pre-registration: [BENCHMARK.md](BENCHMARK.md))
- **Current status and resume point:** [STATUS.md](STATUS.md)
- **Decision log:** [DECISIONS.md](DECISIONS.md)

## Running
Everything runs from one command. It installs Java (if missing) and Nextflow, fetches what each workflow needs and
runs it in containers:

```bash
git clone https://github.com/manuelgug/16s-function-inference && cd 16s-function-inference
./run.sh bench_test   # benchmark smoke test on 1,500 genomes of PICRUSt2's reference (~1 h)
./run.sh bench        # the benchmark: all 27,870 genomes; BENCH_PHASE=dev, later confirm
./run.sh test         # GTDB genome pipeline (annotation, Aim 2), 24 genomes
./run.sh pilot        # GTDB genome pipeline, 2,000 genomes (config/pilot_genomes.tsv)
./run.sh full         # GTDB genome pipeline, every genome passing the filters, or FULL_N of them
```

- **Parameters and primers:** edit the `USER PARAMETERS` block at the top of [`run.sh`](run.sh), or override one for a
  single run, e.g. `./run.sh bench --bench_arms full`.
- **Resuming:** run the same command again; finished steps are skipped.
- **Containers:** Apptainer by default; `CONTAINER=docker ./run.sh bench_test` uses Docker.

### The predictor (s16fun)
```bash
# reference: one 16S per genome, KO copy numbers, GTDB taxonomy (the benchmark builds it from PICRUSt2's reference)
bin/s16fun.py build-ref --ssu ssu.fna --genomes genomes.tsv --ko-copies ko_copies.npz --kos kos.txt \
    --model stacker_all.txt --out ref/
# ASVs (+ optional abundance table) -> per-ASV KO probabilities and calls, copies, modules, novelty,
# and community KO abundances (16S-copy normalised)
bin/s16fun.py predict --ref ref/ --seqs asvs.fna --table counts.tsv --out out/
```
Run inside the analysis image (vsearch, LightGBM): `apptainer exec .apptainer/ghcr.io-manuelgug-16s-function-inference-analysis-0.3.0.img ...`.

### Requirements
- Linux with `curl`, `gzip`, `tar`, `md5sum`; Apptainer or Docker.
- Benchmark: 128 GB RAM (PICRUSt2's hidden-state prediction peaks near 120 GB per fold), ~50 GB disk.
- GTDB pipeline: ~80 GB of databases; 32 GB RAM minimum (eggNOG).

### Outputs
- `results/<profile>/`: everything a run produces, including large tables.
- `reports/<profile>/`: the reports, committed. Each is Markdown plus a self-contained HTML twin (figures embedded and
  explained); open the `.html` in a browser. Benchmark: `reports/bench/dev/benchmark_report.html` (then `confirm/`).
  To render a hand-written report: `apptainer exec .apptainer/ghcr.io-manuelgug-16s-function-inference-analysis-0.3.0.img python bin/render_html.py report.md`.

## Layout
| Path | What |
|---|---|
| `run.sh` | Entry point and all user parameters |
| `benchmark.nf`, `modules/bench.nf` | Benchmark workflow |
| `main.nf`, `modules/` | GTDB genome pipeline (annotation, Aim 2); one file per step |
| `bin/` | Python scripts the steps call (`s16fun.py`: the predictor; `bench_*.py`: the benchmark) |
| `scripts/` | Database download, cost report |
| `conf/` | Containers and resources per step |
| `config/` | Committed inputs: genome lists, *E. coli* reference, checksums, frozen benchmark settings |
| `containers/analysis/` | Analysis image definition, published to GHCR |
| `tests/` | pytest unit tests for `bin/` |

## License
MIT

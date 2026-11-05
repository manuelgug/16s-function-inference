# 16s-function-inference

This project asks whether variation in the 16S rRNA gene can be used to infer the functional repertoire (KEGG KOs and modules, COG categories) of prokaryotic genomes. Mutations are analysed at aligned 16S positions, on full-length genes and on in-silico V4 and V3-V4 amplicons. Reference genomes are GTDB R232 representatives.

- **Design and plan:** [DESIGN.md](DESIGN.md)
- **Current status and resume point:** [STATUS.md](STATUS.md)
- **Decision log:** [DECISIONS.md](DECISIONS.md)

## Running
Everything runs from one command. It installs Java (if missing) and Nextflow, downloads and verifies all databases, and runs the pipeline in containers:

```bash
git clone https://github.com/manuelgug/16s-function-inference && cd 16s-function-inference
./run.sh test     # smoke test, 24 genomes
./run.sh pilot    # 2,000 genomes (fixed list in config/pilot_genomes.tsv)
./run.sh full     # every genome passing the filters, or FULL_N of them
```

- **Parameters and primers:** edit the `USER PARAMETERS` block at the top of [`run.sh`](run.sh), or override one for a single run, e.g. `./run.sh pilot --min_ssu_len 1300`.
- **Resuming:** after an interruption, run the same command again. Downloads continue and finished pipeline steps are skipped.
- **Containers:** Apptainer by default; `CONTAINER=docker ./run.sh test` uses Docker.

### Requirements
- Linux with `curl`, `gzip`, `tar`, `md5sum`.
- Apptainer or Docker.
- About 80 GB of free disk for the databases plus the pilot, more at full scale.
- 32 GB RAM minimum; the eggNOG step is the heaviest.

### Outputs
- `results/<profile>/`: targets (KO, KEGG module, COG matrices), 16S alignment and allele features, QC, and Nextflow reports.
- `reports/<profile>/`: the small summary files, committed to the repo.

## Layout
| Path | What |
|---|---|
| `run.sh` | Entry point and all user parameters |
| `main.nf`, `modules/` | Nextflow workflow; one file per step |
| `bin/` | Python scripts the steps call, run inside the analysis container |
| `scripts/setup_db.sh` | Database download and verification |
| `conf/` | Containers and resources per step |
| `config/` | Committed inputs: genome lists, *E. coli* reference, database checksums |
| `containers/analysis/` | Analysis image definition, published to GHCR |
| `tests/` | pytest unit tests for `bin/` |

## License
MIT

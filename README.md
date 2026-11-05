# 16s-function-inference

This project asks whether variation in the 16S rRNA gene can be used to infer the functional repertoire (KEGG KOs and modules, COG categories) of prokaryotic genomes. Mutations are analysed at aligned 16S positions, on full-length genes and on in-silico V4 and V3-V4 amplicons. Reference genomes are GTDB R232 representatives.

- **Design and plan:** [DESIGN.md](DESIGN.md)
- **Current status and resume point:** [STATUS.md](STATUS.md)
- **Decision log:** [DECISIONS.md](DECISIONS.md)

## Running
Status: not yet implemented (phase P0). Once built, the whole project will run with one command:

```bash
./run.sh test    # smoke test, 24 genomes
./run.sh pilot   # ~2,000 genomes
./run.sh full    # scale decided after the pilot
```

Requirements: Linux, Java ≥ 17, Apptainer (default) or Docker (`CONTAINER=docker ./run.sh ...`), about 100 GB of free disk for the pilot. Databases download automatically.

## License
MIT

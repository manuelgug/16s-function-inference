# Decision log

Every project decision is recorded here, append-only, newest at the bottom. A decision that is later reversed gets a new entry that references the old ID; the old entry is not edited.

**Who**:
- `user`: decided by the user.
- `claude`: decided by Claude. These are open to review; say so and they get a new entry.

| ID | Date | Who | Decision | Alternatives considered | Rationale |
|---|---|---|---|---|---|
| D001 | 2026-11-05 | user | Pursue both aims; prediction (Aim 1) is primary, phylogeny-corrected association (Aim 2) secondary | Prediction only; association only | — |
| D002 | 2026-11-05 | user | 16S input: full-length plus in-silico V4 and V3-V4 amplicons | Full-length only; amplicons only | — |
| D003 | 2026-11-05 | user | Functional targets: KEGG KOs + modules, eggNOG/COG categories | Pfam; curated traits | — |
| D004 | 2026-11-05 | user | Genomes: GTDB species representatives, filtered for full-length 16S | RefSeq complete; AnnoTree precomputed | — |
| D005 | 2026-11-05 | user | Features: aligned 16S positions only | k-mers; embeddings; intragenomic copy variation | — |
| D006 | 2026-11-05 | user | Pilot: ~2,000 genomes stratified across the tree | Single clade; 100 + 5k | — |
| D007 | 2026-11-05 | user | Workflow: Nextflow + containers | Snakemake + conda; Snakemake + Apptainer | — |
| D008 | 2026-11-05 | user | Default container runtime: Apptainer; Docker as a profile | Docker default | — |
| D009 | 2026-11-05 | user | Annotation: eggNOG-mapper only, no KofamScan | eggNOG + KofamScan pilot check; both everywhere | — |
| D010 | 2026-11-05 | user | Analysis language: Python | R; mixed | — |
| D011 | 2026-11-05 | user | Full-scale size decided after the pilot | All filtered reps; capped 20k | — |
| D012 | 2026-11-05 | user | Repo `manuelgug/16s-function-inference`, public, MIT | `16s-mutations`; GPL-3.0 | — |
| D013 | 2026-11-05 | claude | GTDB release R232 (released 2026-04-15) | R226 | Latest release |
| D014 | 2026-11-05 | claude, approved by user | CheckM2 ≥ 90 % complete, ≤ 5 % contaminated; 16S ≥ 1,400 bp. Leaves 54,963 genomes (52,888 bac + 2,075 arc) | Looser thresholds | False gene absences must be rare; full-length 16S is required for the alignment |
| D015 | 2026-11-05 | claude, approved by user | Primers: V4 515F-Y/806RB; V3-V4 341F/805R. Editable in `run.sh` | Other primer pairs | Most common amplicon pairs in current studies |
| D016 | 2026-11-05 | claude, approved by user | CV holdout levels: random, genus, family, order, class | Fewer levels | Map accuracy against phylogenetic distance |
| D017 | 2026-11-05 | claude, approved by user | Filters: KO prevalence 5–95 %; allele frequency ≥ 1 %; alignment column coverage ≥ 95 % | — | Drop trivial targets and noise features |
| D018 | 2026-11-05 | claude, approved by user | Models: null, k-NN, random forest, LightGBM; PICRUSt2 out of the box as a flagged reference | Neural nets; L1 logistic | Cover the baseline (transfer) vs. learned contrast; CPU-only machine |
| D019 | 2026-11-05 | claude, approved by user | Aim 2 covers bacteria only in the pilot | Include archaea | The pilot has only about 80 archaea |
| D020 | 2026-11-05 | claude, approved by user | Analysis container is public on GHCR, built by GitHub Actions | Local build only | Any PC can pull it; no local Docker build needed |
| D021 | 2026-11-05 | user | All scientific parameters and primer sequences editable in one block at the top of `run.sh`; CLI overrides also allowed | Config files only | User request |
| D022 | 2026-11-05 | claude | Git identity in this repo is `manuelgug` with GitHub noreply email (repo-local config) | Work email | No global git identity existed; keeps the work email out of public history |
| D023 | 2026-11-05 | claude | SILVA SSU Ref NR99 **release 144**, from `ftp.arb-silva.de/release_144/` (md5 `fb40d953c6bfef5929a8b7f40fa90184`) | 138.2 | Latest release; the versioned path is stable |
| D024 | 2026-11-05 | claude | Do not spoof a browser user agent to reach SILVA's website | — | Bot-check circumvention is not acceptable; the FTP host works with plain curl |
| D025 | 2026-11-05 | claude | Download genomes from the NCBI FTP tree over HTTPS (`ftp.ncbi.nlm.nih.gov/genomes/all/...`) | `datasets` CLI; Datasets v2 REST API | The REST API returned HTTP 500 on every call (2026-11-05). The newest `datasets` BioContainer is 14.26.0 vs 18.40.0 on conda-forge. The FTP layout is long-stable and needs no extra tool |
| D026 | 2026-11-05 | claude | barrnap **0.9** (`quay.io/biocontainers/barrnap:0.9--hdfd78af_4`) | barrnap 1.10.6 (Infernal-based rewrite) | Both found 7/7 16S copies in *E. coli* K-12. 0.9 is 4× faster (19 s vs 79 s) and is the Prokka/GTDB R232 method |
| D027 | 2026-11-05 | claude | pyseer runs from its own BioContainer (`pyseer:1.4.2--pyhdfd78af_0`), not the analysis image | Analysis image | pyseer 1.4.2 needs `glmnet_py`, which conda cannot solve on linux-64 |
| D028 | 2026-11-05 | claude | Analysis image: `mambaorg/micromamba:2.10.0-ubuntu24.04`, with every package pinned exactly in `containers/analysis/environment.yml` | Unpinned | Reproducibility |
| D029 | 2026-11-05 | claude | *E. coli* reference for position numbering: rrsB, NC_000913.3:4166659-4168200 (1,542 bp), committed as `config/ecoli_rrsB_16S.fna` | J01695 download | Same canonical numbering; no extra download at runtime |
| D030 | 2026-11-05 | claude | Database downloads run on the host from `scripts/setup_db.sh` (resumable `curl -C -`, checksums verified). Only the SINA index build runs as a Nextflow process | Every download as a Nextflow process | A Nextflow task restarts a 12 GB download from zero after interruption. The host only needs curl, gzip, tar and md5sum. Amends DESIGN §7 "DB sub-workflow" |
| D031 | 2026-11-05 | claude | KEGG module completeness from `give_completeness -i <genome\tKOs> -m`, one call for all genomes | Per-genome calls | One process; per-row output is per genome |
| D032 | 2026-11-05 | claude | The decision log is named `DECISIONS.md`, uppercase like `STATUS.md` and `DESIGN.md` | `decisions.md` | Consistent naming. Rename on request |
| D033 | 2026-11-05 | claude | Pilot list: 2,000 genomes, 148 phyla, seed 42, floor 10 per phylum → 1,796 bacteria + 204 archaea; 10 genomes use translation table 4. Test list: 24 genomes from the pilot list, ≥ 2 archaea (got 4), 18 phyla | — | The floor gave more archaea than the ~80 estimated in D019. D019 (Aim 2 bacteria-only in the pilot) stands for now; revisit at P3 |
| D034 | 2026-11-05 | claude | SINA index: a Nextflow `storeDir` process copies the ARB file and lets SINA build `silva.sidx` (458 MB, ~3 min) into `db/silva_<rel>/sina_index/` | Build in `setup_db.sh` | Needs the SINA container; storeDir makes it run once per machine |
| D035 | 2026-11-05 | claude | Added `FULL_N` (0 = all filtered genomes). Pilot and test lists are fixed and committed; the genome filters in `run.sh` affect the full scale only. `MIN_SSU_LEN` also filters our own barrnap copies in every profile | Regenerate pilot list each run | Identical pilot on every PC (DESIGN §2.2) |
| D036 | 2026-11-05 | claude | Nextflow is launched from `results/<profile>/` with work dir `work/<profile>/` | One shared launch dir | `-resume` resumes the last session of the launch dir; per-profile dirs keep test and pilot caches separate |
| D037 | 2026-11-05 | claude | No `subworkflows/` folder; `main.nf` wires the nine modules directly | Subworkflows per stage | The workflow is ~60 lines; another layer adds navigation cost. Amends DESIGN §7 layout |
| D038 | 2026-11-05 | claude | Java: system Java ≥ 17 if present, else `run.sh` downloads a Temurin 21 JRE into `.bin/jre` (Adoptium API) | Require user install | One-script setup on any PC |
| D039 | 2026-11-05 | claude | Nextflow 26.10.4, pinned, self-contained `-dist` jar in `.bin/` | 26.10.5 (3 days old) | Installed and working on this machine; avoid brand-new releases |
| D040 | 2026-11-05 | claude | eggNOG-mapper per batch of 100 genomes, one task at a time with all CPUs, no `--dbmem`, 48 GB task memory | Parallel tasks with `--dbmem` | `--dbmem` loads ~45 GB per task; disk-backed SQLite is enough to start. The pilot timing decides whether to change it |
| D041 | 2026-11-05 | claude | Gene calling: pyrodigal single-genome mode with the GTDB translation table; metagenomic mode if the genome is < 20 kb | Prodigal binary | Same algorithm, Python API, no extra container |
| D042 | 2026-11-05 | claude | NCBI politeness: at most 2 fetch tasks × 4 connections; User-Agent names the project URL, no personal email | More parallelism | Avoid throttling; no personal data in requests |
| D043 | 2026-11-05 | claude | QC "usable" genome = in the 16S feature matrix and ≥ 1 KO annotation; the pilot yield check (≥ 90 %) uses this | Stricter definitions | Minimum needed to enter modelling |
| D044 | 2026-11-05 | claude | Analysis image tag comes from `containers/analysis/VERSION` (0.1.0); bump it whenever `environment.yml` changes | `latest` tag | Pinned, reproducible tag in `conf/base.config` |
| D045 | 2026-11-05 | claude | SINA alignment identity and quality (`align_ident_slv`, `align_quality_slv`) are recorded per genome for QC; no filter on them yet | Filter now | Look at the pilot distribution first |

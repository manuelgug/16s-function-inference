# 16S → function inference: experimental design and plan

**Status:** APPROVED 2026-11-05 (all §10 defaults accepted). Parameters and primer sequences are user-editable at the top of `run.sh`.
**Repo:** https://github.com/manuelgug/16s-function-inference
**Decisions:** every decision, with alternatives and rationale, is logged in [DECISIONS.md](DECISIONS.md). Where this document and the log disagree, the log is newer.

---

## 1. Question

Can the sequence of the 16S rRNA gene, at the level of individual aligned positions ("mutations" / mutation profiles), be used to infer the functional repertoire of a prokaryotic genome?

### Aim 1 (primary): prediction
How accurately can KEGG KO presence/absence, KEGG module completeness and COG functional-category profiles be predicted from 16S sequence?
- How does accuracy decay as the query genome gets phylogenetically farther from the nearest reference genome?
- Does a model over aligned 16S positions beat simple nearest-neighbour transfer (the PICRUSt2-style approach)?
- How much of the signal survives when only a V4 or V3-V4 amplicon is available instead of the full-length gene?

### Aim 2 (secondary): association beyond phylogeny
Which 16S positions or alleles associate with specific functions after controlling for phylogeny, GWAS-style?

### Hypotheses
- **H1.** Accuracy is high when close relatives are in the reference set and decays with phylogenetic distance. Position-based models roughly match 1-NN at short distances. They may degrade more gracefully at long distances. This is testable, not assumed.
- **H2.** Amplicons keep most of the signal of full-length 16S. The loss grows at deeper holdout levels.
- **H3.** Phylogeny explains most of the 16S–function covariation. After correction, residual associations are few and concentrated in positive-control situations, e.g. GC-rich stem alleles vs. hyperthermophily markers. A null result for Aim 2 is an acceptable, informative outcome.

### Central confound
Genomes are not independent: 16S and gene content are both inherited along the same tree. Every evaluation in this design is phylogeny-aware: clade-blocked cross-validation, distance-to-reference curves, and a kinship-corrected LMM for association.

---

## 2. Data

### 2.1 Genome universe
GTDB **R232** (released 2026-04-15): 189,801 bacterial and 10,122 archaeal species representatives.

Filters, applied in this order:
1. `gtdb_representative == t`.
2. CheckM2 completeness ≥ 90 % and contamination ≤ 5 %. Absence of a gene must be meaningful.
3. GTDB `ssu_length` ≥ 1,400 bp. This is a metadata prefilter for a full-length 16S.
4. After our own extraction (§3), at least one non-partial 16S copy ≥ 1,400 bp.

Counts after filters 1–3 (R232 metadata, checked 2026-11-05): **52,888 bacteria + 2,075 archaea = 54,963 genomes**. This is the ceiling for the full-scale run.

Genome sequences are fetched from the NCBI FTP tree over HTTPS (`ftp.ncbi.nlm.nih.gov/genomes/all/...`) by assembly accession, with the GTDB `RS_`/`GB_` prefix removed (D025). This avoids the 192 GB / 132 GB GTDB tarballs for the pilot, and the same code path scales.

### 2.2 Pilot sample (~2,000 genomes)
- Stratified across GTDB phyla. Allocation is proportional to the number of filtered genomes per phylum, with a floor of `min(n_phylum, 10)`. Genomes within a phylum are drawn at random with a fixed seed.
- The selection script runs once. Its output, `config/pilot_genomes.tsv`, is committed, so every PC runs the identical pilot.
- Smoke test (`-profile test`): 24 genomes from the pilot list. It checks plumbing only. It still needs the full databases.

### 2.3 Reference databases
Downloaded by `run.sh`. Versions are pinned. Checksums are recorded in `db/MANIFEST.tsv`.

| Resource | Version | Use | Approx. size |
|---|---|---|---|
| GTDB metadata, taxonomy, `bac120`/`ar53` trees | R232 | sampling, taxonomy, phylogeny | ~0.5 GB |
| eggNOG-mapper data (eggNOG 5.0) | emapper 2.1.15 | KO + COG annotation | ~50 GB (verify at setup) |
| SILVA SSU Ref NR99 ARB file | release 144 (D023) | SINA reference alignment | ~1 GB |
| KEGG module definitions | bundled in `kegg-pathways-completeness` 1.4.4 | module completeness | — |
| PICRUSt2 default reference | bundled in PICRUSt2 2.6.4 | external baseline only | ~1 GB |
| NCBI genome assemblies | per accession | input | ~4 MB/genome |

No data is downloaded from KEGG directly, for licence reasons. SILVA is CC-BY 4.0.

---

## 3. Per-genome processing

1. **Fetch** genome FASTA (NCBI FTP over HTTPS, D025).
2. **Gene calling**: `pyrodigal` in single-genome mode, with the GTDB-reported translation table (e.g. code 4 for Mycoplasmatota).
3. **16S extraction**: `barrnap` 0.9 (D026) with `--kingdom bac|arc` from the GTDB domain. Keep non-partial copies ≥ 1,400 bp. Record copy number as metadata. The **representative copy** is the medoid of the genome's copies (minimum summed pairwise distance), with ties going to the longest copy. Copy-variation features are out of scope (decision 2026-11-05).
4. **Annotation**: `eggNOG-mapper` 2.1.15 in DIAMOND mode with default sensitivity. Columns used: `KEGG_ko`, `COG_category`.
5. **Targets** (genome × target matrices):
   - **KO presence/absence.** Only KOs with prevalence between 5 % and 95 % among the analysed genomes are modelled.
   - **KEGG module completeness** (0–1) from each genome's KO set via `kegg-pathways-completeness`. Modules with near-zero variance are dropped.
   - **COG category profile**: fraction of annotated proteins per COG letter. Multi-letter assignments are split equally.
6. **Alignment**: `SINA` against SILVA Ref NR99. This gives one unified bacteria + archaea column space. Columns are mapped to *E. coli* K-12 (J01695) numbering. Columns with < 95 % non-gap coverage among genomes are dropped.
7. **Features: aligned positions**, the representation chosen 2026-11-05.
   - The state per column is one of {A, C, G, U, gap}. IUPAC ambiguity codes count as missing.
   - Each feature is binary: "position *p* carries allele *X*". A "mutation" is a non-consensus allele at *p*.
   - Alleles with frequency < 1 % and invariant columns are dropped.
   - **Amplicon arms**: in-silico PCR with `cutadapt` linked primers (error rate 0.1) on the representative copy:
     - V4 with 515F-Y / 806RB.
     - V3-V4 with 341F / 805R.

     Amplicon features are the full-length features restricted to the alignment columns the amplicon covers. Genomes failing PCR are excluded from that arm only. Primer coverage per clade is reported as a by-product.
8. **Phylogeny**: GTDB `bac120`/`ar53` trees, pruned to the analysed genomes. It is used for:
   - cross-validation groups (via GTDB taxonomy);
   - distance from each test genome to the nearest training genome;
   - the pyseer kinship matrix.

---

## 4. Aim 1: prediction (superseded 2026-11-06 by §B)

The pilot (P3) tested Aim 1 on GTDB R232 genomes with SINA-aligned allele features, k-NN, random forest and
LightGBM (DECISIONS D046–D076; results in git history, `reports/pilot/pilot_report.md`). Its PICRUSt2 comparison
was not fair: different reference sets and sizes, genome-level bootstrap, untuned settings. §B replaces it with a
benchmark on PICRUSt2's own reference.

---

## B. Benchmark against PICRUSt2 (from 2026-11-06)

**Goal (user, 2026-11-06):** know exactly when and why our method beats PICRUSt2, by how much, with a comparison
that is as fair as possible (same genomes, no leakage, broad phylogenetic coverage), and make the method fast.

### B.1 Data: PICRUSt2's own reference
PICRUSt2 2.6.4 ships its reference (PICRUSt2-SC, GTDB r214): 26,868 bacteria + 1,002 archaea, one aligned full-length
16S per genome, the 16S tree, HMM and RAxML model, and KO copy-number tables (eggNOG-mapper annotation). Both
PICRUSt2 and our methods learn from exactly these genomes, sequences and labels. Truth for a held-out genome is its
row of that table (presence = copy number >= 1). KOs scored: present in >= 10 genomes (10,314).
No downloads and no eggNOG are needed (`bin/bench_extract.py`).

### B.2 Splits, lockbox, leakage
- 5-fold clade-blocked CV at random (species), genus, family, order and class (`bin/make_splits.py`); identical-16S
  genomes share a fold; taxon folds are size-balanced.
- **Lockbox:** families are split once, by phylum, into **dev** (method design, tuning) and **confirm** (scored once,
  frozen settings). Training sets always contain every genome outside the held-out fold.
- **Every method gets the same task:** `bin/bench_prepare.py` writes the training list and one query FASTA per
  level x fold x arm; leak checks (query in training, shared taxon) fail the run.
- **PICRUSt2 per fold:** custom reference built from its own files, restricted to the training genomes (alignment
  rows, pruned tree, trait tables; HMM rebuilt from the subset; constant-0 traits dropped because castor's
  max-parsimony fails on them). Placement once per task, then every HSP setting (mp, emp_prob, pic,
  subtree_average) on the placed trees. Remaining bias in PICRUSt2's favour: tree topology and RAxML model were
  estimated with all genomes.

### B.3 Methods
| method | description |
|---|---|
| **s16fun (ours)** | vsearch top-25 hits -> evidence per (query, KO): 1/5/25-NN averages, KO prevalence in the nearest genome's genus..phylum, global prevalence, KO conservation within genera, hit identities -> LightGBM stacker (trained on dev, cross-fitted by family) -> expected-F1 calls (`bin/s16fun.py`) |
| PICRUSt2 2.6.4 | all 4 HSP settings; thresholds of continuous outputs tuned on dev; comparator = its best setting per level, chosen after seeing the scores |
| k-mer neural net | embeRNA design (1-5-mer frequencies, 2 x 512 ReLU, sigmoid per KO), trained per fold and region |
| nearest-neighbour family | nearest genome; 5 nearest; nearest if >= 97 % (Piphillin / Tax4Fun2-like); genus core if >= 97 % (MicFunPred / PanFP-like) |
| null | KOs in >= 50 % of reference genomes |
Not run as released tools: Tax4Fun2 (no maintained distribution), Piphillin (web server only), embeRNA (EC labels
only), MicFunPred (own IMG reference); each is represented by its method on the shared reference.

### B.4 Endpoints and statistics
Pre-registered in [BENCHMARK.md](BENCHMARK.md) before the confirm run: primary endpoint, win rule, secondary
endpoints, cluster bootstrap over genera (genomes of a genus are not independent), Holm correction.
"When and why" analyses: accuracy by novelty (16S identity, NSTI), by phylum, per KO (AUROC), amplicon arms,
ordination with PERMANOVA / ANOSIM / PERMDISP, stacker feature importance and ablations, run time and memory.

### B.5 Real-world test and communities
- **Real-world:** GTDB R232 genomes not in PICRUSt2's reference, annotated with the GTDB pipeline (`main.nf`);
  PICRUSt2 out of the box against s16fun trained on the same 27,870 genomes. Annotation concordance with
  PICRUSt2's tables is checked first.
- **Communities:** simulated from held-out genomes (exact truth); per-sample error against a permuted-label null and
  differential-abundance concordance (Sun et al. 2020).

### B.6 The tool
`bin/s16fun.py build-ref` (16S FASTA + KO copy table + taxonomy) and `predict` (ASV FASTA [+ abundance table] ->
per-ASV KO probabilities and calls, expected copies, module completeness, novelty, community KO abundances). The
benchmark runs the same code.

## 5. Aim 2: association beyond phylogeny (paused 2026-11-06 until the benchmark is done)

- **Tool**: `pyseer` LMM. The kinship matrix comes from the pruned GTDB tree (`phylogeny_distance.py --lmm`).
- **Domains**: bacteria only in the pilot, since the pilot holds about 80 archaea. Archaea are added at full scale if n permits.
- **Genotypes**: the binary allele-at-position features from §3.7, given as a pyseer `--pres` matrix.
- **Phenotypes**: KO presence (binary, prevalence 5–95 %) and module completeness (continuous).
- **Multiple testing**:
  - Bonferroni over the number of unique variant patterns (`count_patterns.py`), per phenotype;
  - plus BH-FDR across all phenotype × variant tests.

  QQ plots and genomic inflation λ are reported per phenotype.
- **Homoplasy**: an allele–function association is only interpretable if the allele arose more than once on the tree. Report the Fitch-parsimony number of state changes per significant allele.
- **Controls**:
  - **Positive**: reverse gyrase (K03170, a hyperthermophily marker) is expected to associate with G/C alleles at paired stem positions. Stem positions come from the *E. coli* secondary structure.
  - **Calibration**: KOs confined to a single clade should show no residual association after the LMM.
- **Pilot role**: Aim 2 in the pilot is a power and calibration check, not a discovery claim.

---

## 6. Phases, pilot and go/no-go

| Phase | Content | Exit criterion |
|---|---|---|
| P0 | Design, repo, status tracking | User approves this document |
| P1 | Scaffold: `run.sh`, Nextflow workflow, containers, DB setup, pilot sample list, smoke test | `./run.sh test` passes from a fresh clone |
| P2 | Pilot data: 2k genomes, then features + targets | Matrices built; QC report |
| P3 | Pilot analysis: Aims 1 and 2 | Figures and tables in `reports/pilot/` |
| P4 | Pilot report and scale decision | User decides scope of full run |
| P5 | Full scale (≤ 54,963 genomes, scope per P4) | Same outputs at scale |
| P6 | Final report | `reports/final/` |
| B0 | Benchmark infrastructure: PICRUSt2 per-fold references, s16fun, baselines, `./run.sh bench_test` | bench_test passes end to end |
| B1 | Dev runs: method development, speed; dev report | User reviews the dev report |
| B2 | Freeze: BENCHMARK.md + frozen stackers / thresholds committed | User approves the frozen configuration |
| B3 | Confirm run (Track A) + simulated communities | Confirm report |
| B4 | Real-world test (GTDB R232 genomes new since r214) | Report section |
| B5 | Final benchmark report; then the P5 scale decision | User decision |

The scale-up decision is deliberately deferred to after the pilot (decision 2026-11-05). The pilot report gives the user:
- **Engineering**: whether `./run.sh pilot` ran end to end from a clean clone with no manual steps, and whether `-resume` worked after an interruption.
- **Cost**: measured wall time, CPU-hours and disk per genome for each step, extrapolated to 20k and 55k genomes.
- **Yield**: the fraction of sampled genomes with a full-length 16S and an eggNOG annotation. ≥ 90 % is expected; less flags a problem.
- **Signal**: the distance curves, whether M3/M4 beat M1 at any holdout level, the amplicon loss, and Aim 2 calibration (λ and the positive control).

---

## 7. Reproducibility and engineering

- **One entry point**: `./run.sh {test|pilot|full}`.
  - Checks Java ≥ 17 and a container runtime: Apptainer by default, Docker via `CONTAINER=docker`.
  - Installs a pinned Nextflow into `.bin/` if missing.
  - Runs `nextflow run main.nf -profile <scale>,<runtime> -resume`.
  - Database downloads run on the host from `scripts/setup_db.sh` (resumable, checksum-verified, idempotent; D030). The SINA index is built once by a Nextflow process into `db/`.
- **Containers**:
  - Every tool uses a version-pinned BioContainers image.
  - The Python analysis stack (pandas, scikit-learn, LightGBM, pyseer, matplotlib) is one image. It is built from `containers/analysis/Dockerfile` by GitHub Actions, published to GHCR, and referenced by digest.
- **Determinism**: fixed seeds; committed sample list and CV split files.
- **Outputs**:
  - `results/` (large, gitignored) holds per-genome and matrix outputs.
  - `reports/<phase>/` (small, committed) holds tables, figures and the run's `nextflow` report and trace, so results travel with the repo.
- **Tests**: `pytest` on the `bin/` scripts with small fixtures, run in CI. CI does not run the full pipeline.
- **Resources**: `conf/base.config` sizes processes for a 32-core / 125 GB machine. All values can be overridden (`--max_cpus`, `--max_memory`).

### Planned layout
```
.
├── README.md            # what this is, how to run
├── DESIGN.md            # this document
├── STATUS.md            # living task tracker (resume point)
├── run.sh               # single entry point
├── main.nf              # Nextflow workflow
├── nextflow.config
├── conf/                # base / test / pilot / full / apptainer / docker
├── modules/             # one .nf per process
├── scripts/setup_db.sh  # database download (host side)
├── bin/                 # Python scripts called by processes (on PATH in Nextflow)
├── config/              # pilot_genomes.tsv, primers.tsv, splits/
├── containers/analysis/ # Dockerfile for the Python analysis image
├── tests/               # pytest + fixtures
├── reports/             # committed small outputs per phase
├── db/        (ignored) # downloaded databases
├── results/   (ignored)
└── work/      (ignored) # Nextflow work dir
```

---

## 8. Risks and limitations

- **eggNOG-only labels**: KOs come from orthology transfer, not HMM thresholds, and are not cross-checked against KofamScan (decision 2026-11-05). Label noise caps measurable accuracy.
- **Absence ≠ absent**: incomplete genomes and annotations produce false absences. The CheckM2 ≥ 90 % filter mitigates this but does not remove it.
- **One genome per species**: GTDB representatives ignore within-species variation. Near-identical 16S across species bounds the achievable resolution; this is quantified in §4.2.
- **Isolate bias**: MAGs dominate GTDB (151,753 of 189,801 bacterial reps are "derived from metagenome"). MAGs often lack a full-length 16S, so the filter enriches for isolates. Sample composition is reported.
- **Intragenomic 16S heterogeneity**: ignored by the medoid choice.
- **PICRUSt2 leakage**: M2 is optimistic and only a reference point.
- **Disk and time**: the eggNOG DB is about 50 GB. Full-scale cost is unknown until the pilot measures it.

---

## 9. Decision log (Q&A 2026-11-05)

| Topic | Decision |
|---|---|
| Goal | Both aims, predictive first |
| 16S input | Full-length + in-silico V4 and V3-V4 amplicons |
| Functions | KEGG KOs + modules; eggNOG/COG categories |
| Genomes | GTDB representatives, filtered for full-length 16S |
| Features | Aligned positions only |
| Pilot | ~2k genomes stratified across the tree |
| Workflow | Nextflow + containers |
| Default container runtime | Apptainer (Docker as profile) |
| Annotation | eggNOG-mapper only |
| Analysis language | Python |
| Full scale | Decide after pilot |
| Repo | `manuelgug/16s-function-inference`, public, MIT |

## 10. Defaults to confirm at approval

These were chosen by Claude, not asked. Approve or change them.

1. **Quality thresholds**: CheckM2 ≥ 90 % complete, ≤ 5 % contaminated; 16S ≥ 1,400 bp.
2. **Primers**: V4 = 515F-Y / 806RB; V3-V4 = 341F / 805R.
3. **Holdout levels**: random, genus, family, order, class.
4. **Target filters**:
   - KO prevalence 5–95 %;
   - minor-allele frequency ≥ 1 %;
   - alignment column coverage ≥ 95 %.
5. **Models**: k-NN, random forest, LightGBM, plus PICRUSt2 out of the box as a flagged reference.
6. **Aim 2 scope**: bacteria only in the pilot.
7. **Analysis image**: public on GHCR.

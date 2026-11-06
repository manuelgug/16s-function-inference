# Pilot report (P4): 2,000 GTDB R232 genomes

Date: 2026-11-06. Run: `./run.sh pilot`, two Nextflow sessions (8 h 46 m launch, plus a 1 h 05 m resume that added
the steps wired in after launch). Machine: 32 cores, 125 GB RAM. Detailed outputs:

- [qc/qc_report.md](qc/qc_report.md) — yield per step
- [aim1_report/aim1_report.md](aim1_report/aim1_report.md) — prediction (Aim 1)
- [aim2_report/aim2_report.md](aim2_report/aim2_report.md) — association (Aim 2)
- [cost.md](cost.md) — CPU and wall time per step

## 1. Engineering

| check | result |
|---|---|
| tasks failed | 0 (100 in run 1, 89 new + 81 cached in the resume) |
| usable genomes | 1,999 / 2,000 (99.95 %; 1 genome without a full-length 16S). Design target ≥ 90 % |
| in-silico PCR success | V4 98.0 %, V3V4 95.9 % (archaea 90.7 % / 77.5 %) |
| disk | `db/` 51 GB (eggNOG 48 GB), `work/pilot` 9.5 GB, `results/pilot` 228 MB |
| peak memory | PICRUSt2 61.5 GB; eggNOG 11.6 GB; ALIGN_16S 10.3 GB; AIM1_PREDICT 7.4 GB; others < 4 GB |

Problems found and fixed during the pilot:

- **BUILD_TARGETS took 2.3 h on one core**, all in `give_completeness` (~4.7 s per genome; ~70 h at 55k).
  Now runs on 16 genome chunks in parallel, output identical to serial (D067).
- The resume re-ran ALIGN_16S and its dependants once, because deep caching was enabled after the pilot was
  launched (the hash mode changed). Not a bug: the test profile resumes with full caching. `cost_report.py` now counts
  each task once (D068).

## 2. Cost extrapolation

Per-genome steps cost **283 CPU-s per genome**, 98 % of it eggNOG-mapper (runs one batch at a time on all cores,
26.5 cores busy on average).

| step | pilot (2k) | 20k genomes | 54,963 genomes |
|---|---|---|---|
| eggNOG (elapsed, serial batches) | 5.8 h | ~58 h | **~160 h (6.7 days)** |
| fetch + barrnap + genes (run alongside eggNOG) | 2.4 CPU-h | ~24 CPU-h | ~66 CPU-h |
| BUILD_TARGETS (after D067) | ~10 min | ~1.5 h | ~4 h |
| Aim 1 (15 arm × level tasks) | 6.6 CPU-h | ~70–150 CPU-h | ~200–400 CPU-h |
| Aim 2 (pyseer LMM, as built) | 8.8 CPU-h | see §5 | **infeasible as built** (see §5) |
| PICRUSt2 | 13 min, 61.5 GB | needs chunking | needs chunking |
| disk, `work/` | 9.5 GB | ~95 GB | ~260 GB |

Aim 1 is extrapolated from its parts: the k-NN identity matrix (test × train, per fold) grows with n², and the forest grows
with about n log n. At 55k the largest identity block is ~11k × 44k (2 GB float32), so it fits.

## 3. Aim 1 — can 16S predict function?

**Yes. The signal is phylogenetic, and it fades with distance to the nearest reference.**

- Gain in per-genome KO-set F1 over the null model (training-set prevalence), random forest, full-length 16S:
  random split +0.23, genus holdout +0.20, family +0.16, order +0.12, class +0.07.
  All 95 % CIs exclude 0.
- Median per-KO AUROC, full-length random forest: 0.91 (random) → 0.85 (family) → 0.70 (class).
- Accuracy rises steadily with 16S identity to the nearest training genome (fig 1). KO F1 is ~0.60 below 85 %
  identity and ~0.90 at ≥ 99 %. Tree distance and 16S identity agree (Spearman −0.79).
- **Learned models do not beat 5-NN at this scale.** RF and LightGBM are worse than 5-NN at the random, genus and family
  levels, by 0.01–0.03 in F1 and MAE. At order and class they are equal or slightly better: COG Bray-Curtis +0.004–0.008,
  KO F1 not significant. The 16S works as a phylogenetic placement. Individual positions do not add usable
  information beyond "who is my closest relative" — at least not with 2k genomes.
- **Amplicons cost 9–35 % of the gain.** RF KO F1 gain is 0.226 / 0.210 / 0.206 (full / V3V4 / V4) at the random level,
  and 0.073 / 0.055 / 0.060 at class. V3V4 and V4 are indistinguishable.
- **PICRUSt2 matches 1-NN when it is fair.** On the 1,264 genomes not in PICRUSt2's own reference, its KO F1 is 0.79
  (median NSTI 0.26). Our 1-NN scores 0.79 and 5-NN 0.80 at the random level. On the 735 genomes in its reference it
  scores 0.94, which is leakage, not skill.
- Note: the null model's median per-KO AUROC is below 0.5 (0.38–0.47). The null predicts one constant per fold; when a
  held-out clade carries a KO, the training prevalence is lower. Pooled over folds, that anti-correlates. The null's
  per-genome scores, used for "gain over null", are unaffected.

## 4. Aim 2 — 16S alleles associated with function beyond phylogeny

1,795 bacteria, 3,252 alleles (3,106 unique patterns), 4,616 phenotypes (4,233 KOs, 383 modules), pyseer LMM with
GTDB-tree kinship.

- **Calibration is good.** Median genomic inflation λ is 0.94 (KOs) and 0.99 (modules), with 10–90 % range 0.6–1.2.
- **The positive control works.** Reverse gyrase (K03170, hyperthermophily) has 187 Bonferroni hits, including the top
  12 overall. Its positively associated alleles are 74 % G/C, against 47 % background: the known GC–thermostability
  relation of rRNA stems, recovered without being told.
- 11,591 Bonferroni hits (2,518 phenotypes, 2,272 alleles); 20,801 at BH-FDR < 0.05. 10.7 % of tests were flagged by
  pyseer and excluded.
- **The homoplasy filter does not discriminate.** All 2,272 hit alleles change more than once on the tree (Fitch).
  At 1,795 genomes almost every allele is homoplasic, so "> 1 change" separates nothing. Many hits are likely still
  clade-tagging that the kinship term absorbs only partly. A stricter criterion is needed before reading individual hits
  as biology (§6, question 3).

## 5. What must change before full scale

1. **Aim 2 cannot run at 52k bacteria as built.** Each phenotype runs pyseer separately. Each run re-reads an
   n × n kinship (≈ 30 GB as text at 52k) and re-decomposes it (O(n³)). Even with pyseer's `--save-lmm` / `--load-lmm`
   (decompose once), each phenotype costs O(n² × variants): about 1 CPU-h per phenotype at 52k, or ~4,600 CPU-h in total.
   Options: run Aim 2 on a stratified subsample (~15k bacteria: ~5 min per phenotype, ~12 h on 32 cores), plus
   `--save-lmm`.
2. **PICRUSt2 used 61.5 GB for 2k queries.** Its placement step should run in chunks (sequences are independent),
   e.g. 2k per task.
3. **Patristic distances** use a dense n × n matrix (D063). At 55k that is 12 GB, and the kinship it is built from is
   larger still. It needs a nearest-relative-on-tree computation instead of all pairs.

None of these touches the per-genome steps, which take most of the time.

## 6. Open decisions for the full run

1. **Scale.** All 54,963 filtered genomes: ~7 days of eggNOG, ~260 GB of work disk. Or a stratified 20k subset:
   ~2.5 days. The case for 55k is that Aim 1 accuracy depends on reference density (fig 1). The full set is also
   the only way to test whether learned models beat k-NN when data is plentiful.
2. **Aim 2 sample** (§5.1): stratified subsample of about 15k bacteria, or another size.
3. **Aim 2 homoplasy criterion** (§4). Options: require ≥ 5 independent origins of the allele, or require the association
   to replicate within ≥ 2 independent orders (test within clades).
4. Still open from the design: D048 (fixed hyperparameters, no inner CV), D060 (GC-enrichment control instead of a
   stem map), and D019 (Aim 2 bacteria only; the pilot has 204 archaea).

# Pilot report (P4): 2,000 GTDB R232 genomes

Date: 2026-11-06. Run: `./run.sh pilot`, two Nextflow sessions (8 h 46 m launch, plus a 1 h 05 m resume that added
the steps wired in after launch). Machine: 32 cores, 125 GB RAM. Detailed outputs:

- [QC report](qc/qc_report.md) — yield per step
- [Aim 2 report](aim2_report/aim2_report.md) — association
- [Cost report](cost.md) — CPU and wall time per step

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
| Aim 2 (pyseer LMM, as built) | 8.8 CPU-h | see §6 | **infeasible as built** (see §6) |
| PICRUSt2 | 13 min, 61.5 GB | needs chunking | needs chunking |
| PICRUSt2 comparison | 25 min | ~4 h (linear, not measured) | ~11 h (linear, not measured) |
| disk, `work/` | 9.5 GB | ~95 GB | ~260 GB |

Aim 1 is extrapolated from its parts: the k-NN identity matrix (test × train, per fold) grows with n², and the forest grows
with about n log n. At 55k the largest identity block is ~11k × 44k (2 GB float32), so it fits.

## 3. Aim 1 — can 16S predict function? (superseded)

Aim 1 and its PICRUSt2 comparison were redone as a benchmark on PICRUSt2's own reference (DESIGN §B,
`reports/bench/`). The pilot versions are in git history (commit 57d6ab1). What the pilot showed: accuracy rises with
16S identity to the nearest reference genome; learned models on allele features did not beat 5-NN; amplicons lose
9–35 % of the gain over the null; the head-to-head with PICRUSt2 was a tie overall, on unequal references.

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
  as biology (§7, question 3).

![Aim 2 figure 1](aim2_report/fig1_lambda.png)

> **Calibration of the association test.** One λ per phenotype; λ = 1 means the phylogenetic correction removes the
> inflation caused by shared ancestry, above 1 means leftover structure (false positives), below 1 means conservative.
>
> **Takeaway:** both distributions centre on 1 (65 % of KOs between 0.8 and 1.2), so the test is not inflated by
> ancestry. The KO tail below 0.5 means a conservative test for those KOs: less power, not false positives.

![Aim 2 figure 2](aim2_report/fig2_hits_along_16s.png)

> **Where the significant alleles sit along the 16S** (*E. coli* numbering; grey bands = variable regions V1-V9).
>
> **Takeaway:** hits spread along the whole gene and are not enriched in variable regions (43 % of hits, against 48 %
> of tested alleles). The highest peaks sit near position 125 (conserved stretch between V1 and V2), at the end of V2
> (~235) and in V5 (~830-860). Positions with very many hits mark clades rather than single functions, which is why the
> homoplasy criterion (§7, question 3) matters.

## 5. PICRUSt2 head-to-head (superseded)

Replaced by the benchmark (DESIGN §B). The pilot comparison used different reference sets for the two methods
(ours ~1,600 genomes per fold, PICRUSt2 27,870), so its "where we win" findings mixed method and reference effects.

## 6. What must change before full scale

1. **Aim 2 cannot run at 52k bacteria as built.** Each phenotype runs pyseer separately. Each run re-reads an
   n × n kinship (≈ 30 GB as text at 52k) and re-decomposes it (O(n³)). Even with pyseer's `--save-lmm` / `--load-lmm`
   (decompose once), each phenotype costs O(n² × variants): about 1 CPU-h per phenotype at 52k, or ~4,600 CPU-h in total.
   Options: run Aim 2 on a stratified subsample (~15k bacteria: ~5 min per phenotype, ~12 h on 32 cores), plus
   `--save-lmm`.
2. **PICRUSt2 used 61.5 GB for 2k queries.** Resolved in the benchmark design: PICRUSt2 runs per fold (DESIGN §B).

None of these touches the per-genome steps, which take most of the time.

## 7. Open decisions for the full run

1. **Scale.** All 54,963 filtered genomes: ~7 days of eggNOG, ~260 GB of work disk. Or a stratified 20k subset:
   ~2.5 days. The case for 55k is that accuracy depends on reference density. Deferred until the benchmark (B5).
2. **Aim 2 sample** (§6.1): stratified subsample of about 15k bacteria, or another size.
3. **Aim 2 homoplasy criterion** (§4). Options: require ≥ 5 independent origins of the allele, or require the association
   to replicate within ≥ 2 independent orders (test within clades).
4. **PICRUSt2 with our reference.** Decided 2026-11-06: done the other way round, both methods on PICRUSt2's own
   reference (DESIGN §B).
5. Still open from the design: D048 (fixed hyperparameters, no inner CV), D060 (GC-enrichment control instead of a
   stem map), and D019 (Aim 2 bacteria only; the pilot has 204 archaea).

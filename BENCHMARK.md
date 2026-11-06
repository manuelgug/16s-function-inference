# Benchmark pre-registration

**Status:** DRAFT (B1, development). Frozen at B2: from then on this file, `config/bench_frozen/` and the code at the
freeze commit define the confirm run; nothing may change after the confirm run starts.

Design: [DESIGN.md §B](DESIGN.md). Decisions: D078–D088.
The confirm run starts from an empty work directory (Nextflow does not re-run cached tasks when `bin/` scripts change).

## Data and splits
- Genomes, 16S sequences and KO copy numbers: PICRUSt2 2.6.4's own reference (PICRUSt2-SC; 26,868 bacteria + 1,002
  archaea, GTDB r214). Truth: KO present when copy number >= 1. KOs scored: present in >= 10 genomes.
- Folds: 5-fold clade-blocked at genus, family, order and class (plus random = species, reported, not part of the
  rule), `make_splits.py`, seed 42. Identical-16S genomes share a fold.
- Lockbox: families split once by phylum into dev (50 %) and confirm (50 %), seed 42. Dev genomes are used for all
  development and tuning; confirm genomes are scored once, in the confirm run.

## Methods (frozen at B2)
- s16fun at the freeze commit, with `config/bench_frozen/stacker_all.txt`.
- PICRUSt2 2.6.4: settings mp, emp_prob, pic, subtree_average; presence = copy number >= 1, and for each setting
  also >= the dev-tuned threshold in `config/bench_frozen/picrust2_thresholds.tsv`. The comparator at each level is
  the best of these 8 variants on the confirm genomes (chosen after seeing them).
- Baselines as in DESIGN §B.3.

## Co-primary endpoints and win rule (user, 2026-11-06, D088)
Confirm genomes, Track A, full-length 16S; s16fun minus PICRUSt2's best setting; 95 % CIs from resampling genera
(genomes of a genus are not independent).
- **(a) Per-KO AUROC, superiority.** For each KO with >= 5 carriers and >= 5 non-carriers among the held-out genomes of a
  level (folds pooled), the AUROC of each method's score (s16fun: probability; PICRUSt2: predicted copy number). Statistic:
  mean over KOs of s16fun minus PICRUSt2, where PICRUSt2's setting is the one with the highest mean AUROC at that level.
  CI: 200 resamples of genera, every per-KO AUROC recomputed with the resampling weights.
  **Met:** statistic >= +0.05 and CI lower bound > 0.
- **(b) Per-genome KO F1, non-inferiority.** Difference per genome, s16fun minus PICRUSt2's best variant by mean F1.
  CI: 2,000 resamples of genera. **Met:** CI lower bound >= -0.01.
- **Win:** (a) and (b) at each of genus, family, order and class.

Why not F1 superiority (the rule first chosen, D078): on dev folds an oracle that copies the 5 reference genomes with
truly most similar KO sets reaches F1 0.856 at genus level, against 0.832 for PICRUSt2, so no method that predicts from
relatives can gain +0.05 F1 there. F1 superiority per level is reported as a secondary endpoint.

## Secondary endpoints (reported whatever their direction; Holm-adjusted)
KO F1 superiority per level, precision, recall, module completeness error, copy-number error, per-KO AUROC, amplicon arms (V4, V3-V4), accuracy by
novelty (16S identity, NSTI) and phylum, ordination agreement (Mantel, Procrustes, PERMANOVA, ANOSIM, PERMDISP),
run time and memory; Track B and simulated-community endpoints.

# Project status

**Living tracker.** Read this first when resuming, on any machine or in any session. Update it whenever a task changes state, then commit and push. Uncommitted status is lost status.

- **Last updated:** 2026-11-05 by Claude Code (Opus 5.5)
- **Current phase:** P2 — pilot data run (in flight); P3 code in progress
- **Blocking:** nothing

Legend: `[ ]` todo · `[~]` in flight · `[x]` done · `[!]` blocked · `[-]` dropped

Every decision goes into [DECISIONS.md](DECISIONS.md) when it is made.

For in-flight tasks, record machine, start date and how to resume. A long run is resumed with `./run.sh <profile>`; Nextflow `-resume` skips finished work.

## Next action
Pilot is running (`./run.sh pilot`). Meanwhile: finish the fresh-clone test, then write the Aim 1 report and Aim 2. Note for the report: null-model COG Spearman is already ~0.9, so report gain over null.

## Tasks

### P0 — Design
- [x] Clarify requirements with user (Q&A, decisions in DESIGN.md §9)
- [x] Check GTDB R232 availability and filter counts (54,963 candidate genomes)
- [x] Draft DESIGN.md
- [x] Create public repo `manuelgug/16s-function-inference` (MIT)
- [x] User approval of DESIGN.md (2026-11-05: all §10 defaults accepted; parameters + primers must be editable in run.sh)

### P1 — Scaffold
- [x] `run.sh`: USER PARAMETERS block (all params + primers), Java check/auto-JRE, pinned Nextflow 26.10.4, per-profile launch dir (D021, D036, D038, D039)
- [x] `nextflow.config` + `conf/base.config`; profiles test/pilot/full + apptainer/docker
- [x] DB setup `scripts/setup_db.sh`: GTDB R232, SILVA 144, eggNOG 5.0.2 (51 GB on disk)
- [x] eggNOG md5s recorded in `config/db_checksums.tsv`
- [x] Pin SILVA release + record checksum (SILVA 144, D023)
- [x] Analysis image: pinned env, published to GHCR by Actions, public (anonymous pull OK)
- [x] Genome selection; committed `config/pilot_genomes.tsv` (2,000) and `config/test_genomes.tsv` (24) (D033)
- [x] Modules: fetch (NCBI FTP), barrnap 0.9, pyrodigal + 16S medoid, eggNOG-mapper, targets, SINA index + align, features, QC
- [x] pytest for selection, targets, features (9 tests); CI workflow (fixed root user in job container)
- [x] Dev run of test profile (Docker)
- [x] `./run.sh test` passes end to end with Apptainer: 24/24 usable, all 15 Aim 1 tasks. eggNOG 6m15s for 24 genomes (15.6 s/genome on 32 cores)
- [x] `./run.sh test` from a fresh clone passes; QC report byte-identical (D061). **P1 done.**

### P2 — Pilot data (~2k genomes)
- [x] User suggestion: eggNOG without DIAMOND `--iterate`. Tested: 2.0× slower here, KOs near-identical (Jaccard 0.975) → kept `--iterate` (D066). Pilot unchanged.
- [~] `./run.sh pilot` running (started 2026-11-05 on this workstation; log `results/pilot/.nextflow.log`). Expect ~9 h of eggNOG (2,000 × 15.6 s). Resume: `./run.sh pilot`
- [ ] Record runtime/disk per step from `results/pilot/pipeline_info/trace_*.tsv`
- [ ] QC report: 16S yield, annotation yield, sample composition

### P3 — Pilot analysis
- [x] CV split script + module (`bin/make_splits.py`), tested
- [x] Aim 1 prediction script + module (`bin/aim1_predict.py`, one task per arm × level), tested; runs inside the pilot
- [x] Aim 1 report (`bin/aim1_report.py` + `AIM1_REPORT`): accuracy vs identity, gain over null by level, paired vs k-NN, amplicon arms; checked on test run
- [x] M2 PICRUSt2 step + evaluation, leak flagged per genome (D064)
- [x] GTDB-tree distance to nearest training genome (Fig 4 in Aim 1 report) (D063)
- [x] Aim 2 code (`bin/aim2_prep.py`, `bin/aim2_report.py`, `modules/aim2.nf`): tested; test profile runs end to end (n=20 → 95 % of tests flagged by pyseer, as expected)
- [ ] Commit `reports/pilot/`

### P4 — Pilot report and scale decision
- [ ] Pilot report (engineering, cost extrapolation, yield, signal)
- [ ] User decides full-scale scope

### P5 / P6 — Full scale and final report
- [ ] Defined after P4

## Log
Append-only. Newest first.

- 2026-11-05: P3 code done and tested on test profile (Aim 1 + report, PICRUSt2 M2, tree distances, Aim 2, cost script). Benchmarked eggNOG `--dmnd_iterate no`: 2× slower, kept default (D066). Pilot eggNOG in progress.
- 2026-11-05: `./run.sh test` passes end to end (Apptainer). Fixed SINA CSV header issue and Nextflow value channels. Aim 1 wired in. Pilot started. Decisions D046–D053 (D048 deviates from the design: fixed hyperparameters in pilot).
- 2026-11-05: Scaffold written and pushed. Test-profile dev run (Docker) works through SINA alignment; eggNOG DB still downloading. Decisions D033–D045.
- 2026-11-05: Design approved. P1 started. Picked SILVA 144, barrnap 0.9, NCBI FTP download, separate pyseer container (DECISIONS D023–D031). Analysis image builds locally. Added DECISIONS.md.
- 2026-11-05: Project started. Q&A done. DESIGN.md drafted. Repo created. Waiting on approval.

# Project status

**Living tracker.** Read this first when resuming, on any machine or in any session. Update it whenever a task changes state, then commit and push. Uncommitted status is lost status.

- **Last updated:** 2026-11-05 by Claude Code (Opus 5.5)
- **Current phase:** P1 — scaffold (in flight)
- **Blocking:** nothing

Legend: `[ ]` todo · `[~]` in flight · `[x]` done · `[!]` blocked · `[-]` dropped

Every decision goes into [DECISIONS.md](DECISIONS.md) when it is made.

For in-flight tasks, record machine, start date and how to resume. A long run is resumed with `./run.sh <profile>`; Nextflow `-resume` skips finished work.

## Next action
Continue P1: write `scripts/setup_db.sh` and start the DB download, then `run.sh` and the Nextflow modules.

## Tasks

### P0 — Design
- [x] Clarify requirements with user (Q&A, decisions in DESIGN.md §9)
- [x] Check GTDB R232 availability and filter counts (54,963 candidate genomes)
- [x] Draft DESIGN.md
- [x] Create public repo `manuelgug/16s-function-inference` (MIT)
- [x] User approval of DESIGN.md (2026-11-05: all §10 defaults accepted; parameters + primers must be editable in run.sh)

### P1 — Scaffold
- [ ] `run.sh`: checks Java and the container runtime, installs pinned Nextflow, dispatches profile
- [ ] `nextflow.config` + `conf/` (base, test, pilot, full, apptainer, docker)
- [~] DB setup: `scripts/setup_db.sh` (GTDB R232, eggNOG data, SILVA 144) + `db/MANIFEST.tsv` (D030)
- [x] Pin SILVA release + record checksum (SILVA 144, D023)
- [~] Analysis container: image builds locally with pinned env (D028); GHCR build workflow still to do
- [ ] Genome selection script; commit `config/pilot_genomes.tsv` (seeded)
- [ ] Per-genome modules: fetch, pyrodigal, barrnap + medoid, eggNOG-mapper
- [ ] Feature modules: SINA alignment, E. coli mapping, allele matrix, in-silico PCR
- [ ] Target modules: KO matrix, module completeness, COG profiles
- [ ] pytest fixtures + CI
- [ ] `./run.sh test` passes from a fresh clone

### P2 — Pilot data (~2k genomes)
- [ ] Run pilot data stage; record runtime/disk per step
- [ ] QC report: 16S yield, annotation yield, sample composition

### P3 — Pilot analysis
- [ ] CV split files (random, genus, family, order, class)
- [ ] Aim 1: M0–M4 × full/V4/V3-V4; distance curves
- [ ] Aim 2: pyseer LMM, controls, homoplasy
- [ ] Commit `reports/pilot/`

### P4 — Pilot report and scale decision
- [ ] Pilot report (engineering, cost extrapolation, yield, signal)
- [ ] User decides full-scale scope

### P5 / P6 — Full scale and final report
- [ ] Defined after P4

## Log
Append-only. Newest first.

- 2026-11-05: Design approved. P1 started. Picked SILVA 144, barrnap 0.9, NCBI FTP download, separate pyseer container (DECISIONS D023–D031). Analysis image builds locally. Added DECISIONS.md.
- 2026-11-05: Project started. Q&A done. DESIGN.md drafted. Repo created. Waiting on approval.

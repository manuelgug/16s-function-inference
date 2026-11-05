# Project status

**Living tracker.** Read this first when resuming, on any machine or in any session. Update it whenever a task changes state, then commit and push. Uncommitted status is lost status.

- **Last updated:** 2026-11-05 by Claude Code (Opus 5.5)
- **Current phase:** P1 — scaffold (in flight)
- **Blocking:** nothing

Legend: `[ ]` todo · `[~]` in flight · `[x]` done · `[!]` blocked · `[-]` dropped

Every decision goes into [DECISIONS.md](DECISIONS.md) when it is made.

For in-flight tasks, record machine, start date and how to resume. A long run is resumed with `./run.sh <profile>`; Nextflow `-resume` skips finished work.

## Next action
Wait for eggNOG DB download, then run `./run.sh test` (Apptainer) end to end; fix what breaks; then fresh-clone test.

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
- [~] DB setup `scripts/setup_db.sh`: GTDB R232 + SILVA 144 done and verified; eggNOG 5.0.2 downloading (started 2026-11-05, this workstation). Resume: rerun `scripts/setup_db.sh db 232 144` or `./run.sh test`
- [ ] Record eggNOG md5s in `config/db_checksums.tsv` once the download finishes (script prints them)
- [x] Pin SILVA release + record checksum (SILVA 144, D023)
- [x] Analysis image: pinned env, published to GHCR by Actions, public (anonymous pull OK)
- [x] Genome selection; committed `config/pilot_genomes.tsv` (2,000) and `config/test_genomes.tsv` (24) (D033)
- [x] Modules: fetch (NCBI FTP), barrnap 0.9, pyrodigal + 16S medoid, eggNOG-mapper, targets, SINA index + align, features, QC
- [x] pytest for selection, targets, features (9 tests); CI workflow (fixed root user in job container)
- [~] Dev run of test profile (Docker): fetch 24/24, 16S 24/24, SINA alignment OK; eggNOG waits for DB
- [ ] `./run.sh test` passes end to end with Apptainer (default runtime)
- [ ] `./run.sh test` passes from a fresh clone (then P1 done)

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

- 2026-11-05: Scaffold written and pushed. Test-profile dev run (Docker) works through SINA alignment; eggNOG DB still downloading. Decisions D033–D045.
- 2026-11-05: Design approved. P1 started. Picked SILVA 144, barrnap 0.9, NCBI FTP download, separate pyseer container (DECISIONS D023–D031). Analysis image builds locally. Added DECISIONS.md.
- 2026-11-05: Project started. Q&A done. DESIGN.md drafted. Repo created. Waiting on approval.

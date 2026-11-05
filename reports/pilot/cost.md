# Cost report

Measured on 2000 genomes. Wall hours are summed task time, not elapsed time (tasks overlap).

| step | tasks | wall_h | cpu_h | scaling | cpu_h_at_20000 | cpu_h_at_54963 |
|---|---|---|---|---|---|---|
| AIM1_PREDICT | 15 | 0.58 | 6.51 | ~n log n (RF) to n^2 (k-NN) | nan | nan |
| ALIGN_16S | 1 | 0.01 | 0.10 | ~linear | nan | nan |
| BARRNAP | 20 | 0.22 | 1.42 | linear (per genome) | 14.15 | 38.89 |
| BUILD_FEATURES | 1 | 0.00 | 0.00 | ~linear | nan | nan |
| BUILD_TARGETS | 1 | 2.30 | 2.30 | ~linear | nan | nan |
| EGGNOG | 20 | 5.83 | 154.57 | linear (per genome) | 1545.68 | 4247.76 |
| FETCH_GENOMES | 20 | 0.63 | 0.01 | linear (per genome) | 0.11 | 0.30 |
| MAKE_SPLITS | 1 | 0.00 | 0.00 | small | nan | nan |
| PROCESS_GENOMES | 20 | 0.16 | 0.98 | linear (per genome) | 9.77 | 26.86 |
| QC_REPORT | 1 | 0.00 | 0.00 | small | nan | nan |

Per-genome steps: 282.5 CPU-seconds per genome.
- 20000 genomes: ~1570 CPU-hours for per-genome steps (~49 h on 32 cores).
- 54963 genomes: ~4314 CPU-hours for per-genome steps (~135 h on 32 cores).

# Cost report

Measured on 2000 genomes. Wall hours are summed task time, not elapsed time (tasks overlap).

| step | tasks | wall_h | cpu_h | scaling | cpu_h_at_20000 | cpu_h_at_54963 |
|---|---|---|---|---|---|---|
| AIM1_PREDICT | 15 | 0.61 | 6.64 | ~n log n (RF) to n^2 (k-NN) | nan | nan |
| AIM1_REPORT | 1 | 0.00 | 0.00 | small | nan | nan |
| AIM2_PREP | 1 | 0.02 | 0.02 | small | nan | nan |
| AIM2_PYSEER | 64 | 8.82 | 8.77 | ~n^2 per test (LMM) | nan | nan |
| AIM2_REPORT | 1 | 0.01 | 0.01 | small | nan | nan |
| ALIGN_16S | 1 | 0.01 | 0.09 | ~linear | nan | nan |
| BARRNAP | 20 | 0.22 | 1.42 | linear (per genome) | 14.15 | 38.89 |
| BUILD_FEATURES | 1 | 0.00 | 0.00 | ~linear | nan | nan |
| BUILD_TARGETS | 1 | 2.30 | 2.30 | ~linear | nan | nan |
| EGGNOG | 20 | 5.83 | 154.57 | linear (per genome) | 1545.68 | 4247.76 |
| FETCH_GENOMES | 20 | 0.63 | 0.01 | linear (per genome) | 0.11 | 0.30 |
| MAKE_SPLITS | 1 | 0.00 | 0.00 | small | nan | nan |
| PATRISTIC | 1 | 0.02 | 0.02 | n^2 memory (D063) | nan | nan |
| PICRUST2 | 1 | 0.22 | 0.92 | ~linear + fixed reference cost | nan | nan |
| PICRUST2_EVAL | 1 | 0.00 | 0.00 | small | nan | nan |
| PROCESS_GENOMES | 20 | 0.16 | 0.98 | linear (per genome) | 9.77 | 26.86 |
| QC_REPORT | 1 | 0.00 | 0.00 | small | nan | nan |

Per-genome steps: 282.5 CPU-seconds per genome.
- 20000 genomes: ~1570 CPU-hours for per-genome steps (~49 h on 32 cores).
- 54963 genomes: ~4314 CPU-hours for per-genome steps (~135 h on 32 cores).

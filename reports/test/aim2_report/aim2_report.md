# Aim 2 report: 16S alleles associated with function beyond phylogeny

Genomes: 20 (bacteria). Variants tested: 2123 (1070 unique patterns; Bonferroni threshold p < 4.67e-05). Phenotypes: 5093 KOs, 383 modules.

## Calibration

![lambda](fig1_lambda.png)

| kind | count | 10% | 50% | 90% |
|---|---|---|---|---|
| KO | 5.09e+03 | 0.324 | 0.724 | 1.45 |
| module | 383 | 0.68 | 0.986 | 1.25 |

## Significant associations

Bonferroni hits (phenotype × allele): **388** across 126 phenotypes and 220 alleles. BH-FDR < 0.05 (unflagged): 826. Tests flagged by pyseer (`notes`, e.g. bad-chisq / high-bse) are never hits: 95.0% of all tests.

Alleles among hits that arose more than once on the tree (Fitch changes > 1): **204 of 220**. Only these can support association beyond phylogeny.

![hits](fig2_hits_along_16s.png)

Top 30 hits:

| phenotype | variant | ecoli_pos | state | af | beta | p | fitch_changes |
|---|---|---|---|---|---|---|---|
| M00128 | 689_21786_A | 689 | A | 0.1 | 0.125 | 0 | 2 |
| M00128 | 698_21917_T | 698 | T | 0.1 | 0.125 | 0 | 2 |
| M00128 | 974_31158_C | 974 | C | 0.1 | 0.125 | 0 | 2 |
| M00883 | 763_22539_C | 763 | C | 0.1 | -0.333 | 0 | 2 |
| M00883 | 465_9992_G | 465 | G | 0.1 | -0.333 | 0 | 2 |
| M00883 | 578_15646_G | 578 | G | 0.1 | -0.333 | 0 | 2 |
| M00993 | 834_26161_C | 834 | C | 0.15 | 0.5 | 0 | 3 |
| M00993 | 852_26802_G | 852 | G | 0.15 | 0.5 | 0 | 3 |
| M00085 | 128_2062_G | 128 | G | 0.8 | -0.333 | 5.07e-144 | 2 |
| M00085 | 233_5272_C | 233 | C | 0.8 | -0.333 | 5.07e-144 | 2 |
| M00128 | 974_31158_A | 974 | A | 0.9 | -0.125 | 2.11e-142 | 2 |
| M00744 | 705_21929_T | 705 | T | 0.9 | -0.667 | 5.98e-140 | 1 |
| M00744 | 719_21995_T | 719 | T | 0.1 | 0.667 | 5.98e-140 | 1 |
| M00744 | 948_29668_T | 948 | T | 0.1 | 0.667 | 5.98e-140 | 1 |
| M00744 | 705_21929_G | 705 | G | 0.1 | 0.667 | 5.98e-140 | 1 |
| M00744 | 1233_40274_A | 1233 | A | 0.1 | 0.667 | 5.98e-140 | 1 |
| M00744 | 719_21995_C | 719 | C | 0.9 | -0.667 | 5.98e-140 | 1 |
| M00157 | 592_15692_C | 592 | C | 0.1 | -0.988 | 5.32e-20 | 2 |
| M00157 | 1024_32538_T | 1024 | T | 0.1 | -0.988 | 5.32e-20 | 2 |
| M00157 | 232_5270_C | 232 | C | 0.1 | -0.988 | 5.32e-20 | 2 |
| M00051 | 489_10262_T | 489 | T | 0.1 | -0.63 | 1.7e-18 | 2 |
| M00051 | 1470_43028_G | 1470 | G | 0.9 | 0.63 | 1.7e-18 | 2 |
| M00149 | 895_27165_G | 895 | G | 0.75 | 0.781 | 4.65e-11 | 4 |
| M00919 | 1172_37683_T | 1172 | T | 0.1 | 0.388 | 1.46e-10 | 2 |
| M00378 | 612_15963_C | 612 | C | 0.85 | -0.516 | 3.3e-10 | 3 |
| M00378 | 628_16292_G | 628 | G | 0.85 | -0.516 | 3.3e-10 | 3 |
| M00810 | 193_4075_G | 193 | G | 0.15 | 0.127 | 4.31e-10 | 3 |
| M00642 | 1233_40274_A | 1233 | A | 0.1 | 0.621 | 4.65e-10 | 1 |
| M00642 | 719_21995_T | 719 | T | 0.1 | 0.621 | 4.65e-10 | 1 |
| M00642 | 705_21929_G | 705 | G | 0.1 | 0.621 | 4.65e-10 | 1 |

## Controls

_no control phenotype present in >= 5 genomes_

Positive control expectation: reverse gyrase (K03170) hits with beta > 0 are G/C-enriched relative to the background allele composition.

Files: `per_phenotype.tsv`, `hits_bonferroni.tsv`.

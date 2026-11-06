# Aim 1 report: predicting function from 16S

Genomes: 24. Arms: V3V4, V4, full. Levels: random, genus, family, order, class.

## 16S identity to the nearest training genome (full length)

How close each test genome's 16S is to its most similar training genome, per holdout level (quantiles over genomes). Stricter holdouts push test genomes further from the training data.

| level | count | 10% | 50% | 90% |
|---|---|---|---|---|
| class | 24.000 | 0.767 | 0.808 | 0.835 |
| family | 24.000 | 0.767 | 0.810 | 0.862 |
| genus | 24.000 | 0.767 | 0.810 | 0.916 |
| order | 24.000 | 0.758 | 0.810 | 0.860 |
| random | 24.000 | 0.767 | 0.813 | 0.916 |

## Accuracy vs 16S identity (full length, all levels pooled)

![fig1](fig1_accuracy_vs_identity.png)

> **Figure 1. Accuracy against 16S similarity to the closest training genome.** Every test genome is placed in a bin by the full-length 16S identity to its most similar genome in the training folds (all holdout levels pooled, so each genome counts once per level; n under each bin). Each point is the mean per-genome score in that bin.
>
> Panels: **KO set F1** = overlap between the predicted and the true set of KEGG orthologs of a genome (1 = perfect); **module completeness MAE** = mean absolute error of predicted KEGG module completeness (0-1 scale; 0 = perfect); **COG profile Bray-Curtis** = dissimilarity between predicted and true COG-category composition (0 = identical).
>
> Models: **null** = training-set prevalence or mean, ignores the 16S (dashed grey); **1-NN / 5-NN** = copy the function profile of the 1 or 5 training genomes with the most similar 16S; **random forest / LightGBM** = learned from aligned 16S positions (LightGBM is fitted for modules and COGs only).
>
> **How to read it:** the distance between a model line and the dashed null line is what the 16S adds. Lines that slope towards perfect (up for F1, down for the errors) mean predictions are only as good as the nearest relative in the reference.
>
> **Here:** the best KO F1 goes from 0.55 in the lowest-identity bin (<85 %) to nan at >=99 %, while the null stays between 0.49 and 0.55.

## Accuracy vs GTDB-tree distance (full length)

Spearman correlation between 16S identity and tree distance to the nearest training genome: -0.789.

![fig4](fig4_accuracy_vs_tree_distance.png)

> **Figure 4. The same accuracy curves against evolutionary distance on the GTDB tree.** The x axis is the patristic distance (sum of branch lengths, substitutions per site on the GTDB marker-protein tree) from each test genome to its closest training genome, cut into five equal-count bins (quintiles). Archaea and bacteria sit on separate trees, so each genome's nearest relative is searched within its domain.
>
> Panels: **KO set F1** = overlap between the predicted and the true set of KEGG orthologs of a genome (1 = perfect); **module completeness MAE** = mean absolute error of predicted KEGG module completeness (0-1 scale; 0 = perfect); **COG profile Bray-Curtis** = dissimilarity between predicted and true COG-category composition (0 = identical).
>
> Models: **null** = training-set prevalence or mean, ignores the 16S (dashed grey); **1-NN / 5-NN** = copy the function profile of the 1 or 5 training genomes with the most similar 16S; **random forest / LightGBM** = learned from aligned 16S positions (LightGBM is fitted for modules and COGs only).
>
> **How to read it:** if this figure looks like figure 1 mirrored, 16S identity is mostly a proxy for phylogenetic distance. Spearman correlation of the two distances: -0.79.
>
> **Here:** the best KO F1 falls from 0.67 (closest quintile) to 0.50 (most distant).

## Gain over the null model by holdout level (full length)

![fig2](fig2_gain_over_null_by_level.png)

> **Figure 2. Gain over the null model, by how distant the test genomes are from the training data.** Each holdout level keeps whole taxa out of training: *random* = ordinary 5-fold split; *genus* ... *class* = every genome of the test genus (or family, order, class) is unseen during training. Points are mean per-genome gain over the null (positive = better than ignoring the 16S); bars are 95 % bootstrap CIs.
>
> Panels: **KO set F1** = overlap between the predicted and the true set of KEGG orthologs of a genome (1 = perfect); **module completeness MAE** = mean absolute error of predicted KEGG module completeness (0-1 scale; 0 = perfect); **COG profile Bray-Curtis** = dissimilarity between predicted and true COG-category composition (0 = identical).
>
> Models: **null** = training-set prevalence or mean, ignores the 16S (dashed grey); **1-NN / 5-NN** = copy the function profile of the 1 or 5 training genomes with the most similar 16S; **random forest / LightGBM** = learned from aligned 16S positions (LightGBM is fitted for modules and COGs only).
>
> **How to read it:** the left-to-right decline shows how fast 16S-based prediction degrades for organisms with no close relative in the reference. A point at 0 means no better than the null.
>
> **Here:** random forest KO F1 gain is 0.037 on a random split and 0.015 when whole classes are held out.

| target_type | metric | model | level | mean_gain | ci_lo | ci_hi | n |
|---|---|---|---|---|---|---|---|
| ko | f1 | M1_knn1 | random | 0.059 | -0.003 | 0.129 | 24 |
| ko | f1 | M1_knn1 | genus | 0.060 | -0.000 | 0.125 | 24 |
| ko | f1 | M1_knn1 | family | 0.049 | -0.005 | 0.107 | 24 |
| ko | f1 | M1_knn1 | order | 0.022 | -0.023 | 0.068 | 24 |
| ko | f1 | M1_knn1 | class | 0.034 | -0.018 | 0.093 | 24 |
| ko | f1 | M1_knn5 | random | 0.075 | 0.041 | 0.111 | 24 |
| ko | f1 | M1_knn5 | genus | 0.094 | 0.052 | 0.136 | 24 |
| ko | f1 | M1_knn5 | family | 0.063 | 0.031 | 0.094 | 24 |
| ko | f1 | M1_knn5 | order | 0.062 | 0.025 | 0.103 | 24 |
| ko | f1 | M1_knn5 | class | 0.065 | 0.031 | 0.100 | 24 |
| ko | f1 | M3_rf | random | 0.037 | 0.018 | 0.057 | 24 |
| ko | f1 | M3_rf | genus | 0.042 | 0.022 | 0.064 | 24 |
| ko | f1 | M3_rf | family | 0.043 | 0.025 | 0.061 | 24 |
| ko | f1 | M3_rf | order | 0.025 | 0.008 | 0.042 | 24 |
| ko | f1 | M3_rf | class | 0.015 | 0.005 | 0.026 | 24 |
| module | mae | M1_knn1 | random | 0.016 | -0.005 | 0.039 | 24 |
| module | mae | M1_knn1 | genus | 0.016 | -0.006 | 0.040 | 24 |
| module | mae | M1_knn1 | family | 0.008 | -0.012 | 0.027 | 24 |
| module | mae | M1_knn1 | order | 0.005 | -0.013 | 0.025 | 24 |
| module | mae | M1_knn1 | class | 0.005 | -0.013 | 0.022 | 24 |
| module | mae | M1_knn5 | random | 0.018 | 0.006 | 0.030 | 24 |
| module | mae | M1_knn5 | genus | 0.023 | 0.009 | 0.035 | 24 |
| module | mae | M1_knn5 | family | 0.019 | 0.005 | 0.033 | 24 |
| module | mae | M1_knn5 | order | 0.017 | 0.005 | 0.029 | 24 |
| module | mae | M1_knn5 | class | 0.018 | 0.006 | 0.029 | 24 |
| module | mae | M3_rf | random | 0.014 | 0.007 | 0.021 | 24 |
| module | mae | M3_rf | genus | 0.018 | 0.009 | 0.026 | 24 |
| module | mae | M3_rf | family | 0.014 | 0.006 | 0.023 | 24 |
| module | mae | M3_rf | order | 0.013 | 0.006 | 0.021 | 24 |
| module | mae | M3_rf | class | 0.013 | 0.006 | 0.020 | 24 |
| module | mae | M4_lgbm | random | -0.001 | -0.002 | 0.000 | 24 |
| module | mae | M4_lgbm | genus | -0.001 | -0.005 | 0.002 | 24 |
| module | mae | M4_lgbm | family | 0.000 | -0.001 | 0.002 | 24 |
| module | mae | M4_lgbm | order | 0.001 | -0.002 | 0.004 | 24 |
| module | mae | M4_lgbm | class | -0.000 | -0.005 | 0.007 | 24 |
| cog | bray_curtis | M1_knn1 | random | -0.030 | -0.056 | -0.009 | 24 |
| cog | bray_curtis | M1_knn1 | genus | -0.032 | -0.056 | -0.011 | 24 |
| cog | bray_curtis | M1_knn1 | family | -0.034 | -0.061 | -0.013 | 24 |
| cog | bray_curtis | M1_knn1 | order | -0.036 | -0.053 | -0.019 | 24 |
| cog | bray_curtis | M1_knn1 | class | -0.028 | -0.041 | -0.014 | 24 |
| cog | bray_curtis | M1_knn5 | random | -0.001 | -0.008 | 0.006 | 24 |
| cog | bray_curtis | M1_knn5 | genus | 0.004 | -0.005 | 0.013 | 24 |
| cog | bray_curtis | M1_knn5 | family | 0.001 | -0.007 | 0.010 | 24 |
| cog | bray_curtis | M1_knn5 | order | -0.004 | -0.012 | 0.004 | 24 |
| cog | bray_curtis | M1_knn5 | class | 0.005 | -0.004 | 0.014 | 24 |
| cog | bray_curtis | M3_rf | random | 0.004 | -0.001 | 0.008 | 24 |
| cog | bray_curtis | M3_rf | genus | 0.007 | 0.001 | 0.013 | 24 |
| cog | bray_curtis | M3_rf | family | 0.006 | 0.001 | 0.011 | 24 |
| cog | bray_curtis | M3_rf | order | 0.004 | -0.003 | 0.010 | 24 |
| cog | bray_curtis | M3_rf | class | 0.007 | 0.002 | 0.013 | 24 |
| cog | bray_curtis | M4_lgbm | random | -0.001 | -0.003 | 0.001 | 24 |
| cog | bray_curtis | M4_lgbm | genus | -0.003 | -0.005 | -0.001 | 24 |
| cog | bray_curtis | M4_lgbm | family | -0.001 | -0.004 | 0.002 | 24 |
| cog | bray_curtis | M4_lgbm | order | 0.001 | -0.001 | 0.004 | 24 |
| cog | bray_curtis | M4_lgbm | class | -0.004 | -0.010 | 0.002 | 24 |

## Learned models vs k-NN (paired per genome; positive = learned model better)

Same genomes, same folds: mean per-genome difference between a learned model and a k-NN model, with 95 % bootstrap CI and Wilcoxon signed-rank p. Positive means the learned model is better (errors are sign-flipped). This tests whether learning from individual 16S positions adds anything over copying the closest relatives.

| arm | level | metric | comparison | n | mean_diff | ci_lo | ci_hi | wilcoxon_p |
|---|---|---|---|---|---|---|---|---|
| V3V4 | class | bray_curtis | M3_rf - M1_knn1 | 24 | 0.015 | -0.005 | 0.035 | 0.023 |
| V3V4 | class | f1 | M3_rf - M1_knn1 | 24 | -0.042 | -0.111 | 0.021 | 0.833 |
| V3V4 | class | mae | M3_rf - M1_knn1 | 24 | -0.005 | -0.030 | 0.019 | 0.603 |
| V3V4 | family | bray_curtis | M3_rf - M1_knn1 | 24 | 0.013 | -0.017 | 0.045 | 0.178 |
| V3V4 | family | f1 | M3_rf - M1_knn1 | 24 | -0.036 | -0.117 | 0.036 | 0.944 |
| V3V4 | family | mae | M3_rf - M1_knn1 | 24 | -0.014 | -0.041 | 0.012 | 0.877 |
| V3V4 | genus | bray_curtis | M3_rf - M1_knn1 | 24 | 0.021 | -0.006 | 0.049 | 0.040 |
| V3V4 | genus | f1 | M3_rf - M1_knn1 | 24 | -0.036 | -0.109 | 0.028 | 0.922 |
| V3V4 | genus | mae | M3_rf - M1_knn1 | 24 | -0.013 | -0.044 | 0.012 | 0.833 |
| V3V4 | order | bray_curtis | M3_rf - M1_knn1 | 24 | 0.027 | 0.001 | 0.053 | 0.021 |
| V3V4 | order | f1 | M3_rf - M1_knn1 | 24 | -0.022 | -0.082 | 0.030 | 0.989 |
| V3V4 | order | mae | M3_rf - M1_knn1 | 24 | -0.007 | -0.032 | 0.015 | 0.790 |
| V3V4 | random | bray_curtis | M3_rf - M1_knn1 | 24 | 0.031 | 0.011 | 0.058 | 0.004 |
| V3V4 | random | f1 | M3_rf - M1_knn1 | 24 | -0.028 | -0.093 | 0.029 | 0.747 |
| V3V4 | random | mae | M3_rf - M1_knn1 | 24 | 0.001 | -0.019 | 0.020 | 0.726 |
| V4 | class | bray_curtis | M3_rf - M1_knn1 | 24 | 0.036 | 0.021 | 0.055 | 0.000 |
| V4 | class | f1 | M3_rf - M1_knn1 | 24 | -0.013 | -0.061 | 0.030 | 0.747 |
| V4 | class | mae | M3_rf - M1_knn1 | 24 | 0.015 | -0.000 | 0.029 | 0.089 |
| V4 | family | bray_curtis | M3_rf - M1_knn1 | 24 | 0.047 | 0.027 | 0.072 | 0.000 |
| V4 | family | f1 | M3_rf - M1_knn1 | 24 | 0.003 | -0.051 | 0.053 | 0.509 |
| V4 | family | mae | M3_rf - M1_knn1 | 24 | 0.017 | 0.001 | 0.032 | 0.074 |
| V4 | genus | bray_curtis | M3_rf - M1_knn1 | 24 | 0.047 | 0.025 | 0.075 | 0.000 |
| V4 | genus | f1 | M3_rf - M1_knn1 | 24 | -0.012 | -0.066 | 0.037 | 0.966 |
| V4 | genus | mae | M3_rf - M1_knn1 | 24 | 0.010 | -0.011 | 0.029 | 0.264 |
| V4 | order | bray_curtis | M3_rf - M1_knn1 | 24 | 0.046 | 0.030 | 0.065 | 0.000 |
| V4 | order | f1 | M3_rf - M1_knn1 | 24 | -0.002 | -0.035 | 0.031 | 0.966 |
| V4 | order | mae | M3_rf - M1_knn1 | 24 | 0.010 | -0.006 | 0.024 | 0.178 |
| V4 | random | bray_curtis | M3_rf - M1_knn1 | 24 | 0.041 | 0.021 | 0.067 | 0.000 |
| V4 | random | f1 | M3_rf - M1_knn1 | 24 | -0.024 | -0.079 | 0.027 | 0.643 |
| V4 | random | mae | M3_rf - M1_knn1 | 24 | 0.003 | -0.019 | 0.022 | 0.790 |
| full | class | bray_curtis | M3_rf - M1_knn1 | 24 | 0.036 | 0.023 | 0.048 | 0.000 |
| full | class | f1 | M3_rf - M1_knn1 | 24 | -0.018 | -0.076 | 0.034 | 0.584 |
| full | class | mae | M3_rf - M1_knn1 | 24 | 0.007 | -0.006 | 0.020 | 0.208 |
| full | family | bray_curtis | M3_rf - M1_knn1 | 24 | 0.040 | 0.019 | 0.065 | 0.001 |
| full | family | f1 | M3_rf - M1_knn1 | 24 | -0.006 | -0.063 | 0.048 | 0.509 |
| full | family | mae | M3_rf - M1_knn1 | 24 | 0.006 | -0.009 | 0.023 | 0.406 |
| full | genus | bray_curtis | M3_rf - M1_knn1 | 24 | 0.039 | 0.021 | 0.062 | 0.001 |
| full | genus | f1 | M3_rf - M1_knn1 | 24 | -0.018 | -0.076 | 0.040 | 0.966 |
| full | genus | mae | M3_rf - M1_knn1 | 24 | 0.002 | -0.015 | 0.019 | 0.663 |
| full | order | bray_curtis | M3_rf - M1_knn1 | 24 | 0.040 | 0.025 | 0.056 | 0.000 |
| full | order | f1 | M3_rf - M1_knn1 | 24 | 0.004 | -0.035 | 0.037 | 0.473 |
| full | order | mae | M3_rf - M1_knn1 | 24 | 0.008 | -0.006 | 0.023 | 0.360 |
| full | random | bray_curtis | M3_rf - M1_knn1 | 24 | 0.034 | 0.012 | 0.056 | 0.004 |
| full | random | f1 | M3_rf - M1_knn1 | 24 | -0.022 | -0.084 | 0.038 | 0.922 |
| full | random | mae | M3_rf - M1_knn1 | 24 | -0.003 | -0.023 | 0.017 | 0.922 |
| V3V4 | class | bray_curtis | M3_rf - M1_knn5 | 24 | 0.011 | 0.004 | 0.018 | 0.015 |
| V3V4 | class | f1 | M3_rf - M1_knn5 | 24 | -0.033 | -0.074 | 0.001 | 0.229 |
| V3V4 | class | mae | M3_rf - M1_knn5 | 24 | 0.002 | -0.007 | 0.011 | 0.944 |
| V3V4 | family | bray_curtis | M3_rf - M1_knn5 | 24 | 0.005 | -0.002 | 0.012 | 0.151 |
| V3V4 | family | f1 | M3_rf - M1_knn5 | 24 | -0.020 | -0.053 | 0.007 | 0.768 |
| V3V4 | family | mae | M3_rf - M1_knn5 | 24 | -0.001 | -0.009 | 0.006 | 0.422 |
| V3V4 | genus | bray_curtis | M3_rf - M1_knn5 | 24 | 0.011 | 0.004 | 0.018 | 0.004 |
| V3V4 | genus | f1 | M3_rf - M1_knn5 | 24 | -0.037 | -0.078 | -0.001 | 0.218 |
| V3V4 | genus | mae | M3_rf - M1_knn5 | 24 | -0.003 | -0.012 | 0.007 | 0.663 |
| V3V4 | order | bray_curtis | M3_rf - M1_knn5 | 24 | 0.011 | 0.004 | 0.017 | 0.005 |
| V3V4 | order | f1 | M3_rf - M1_knn5 | 24 | -0.024 | -0.060 | 0.008 | 0.603 |
| V3V4 | order | mae | M3_rf - M1_knn5 | 24 | -0.001 | -0.009 | 0.008 | 0.877 |
| V3V4 | random | bray_curtis | M3_rf - M1_knn5 | 24 | 0.008 | 0.001 | 0.015 | 0.021 |
| V3V4 | random | f1 | M3_rf - M1_knn5 | 24 | -0.021 | -0.043 | 0.000 | 0.240 |
| V3V4 | random | mae | M3_rf - M1_knn5 | 24 | -0.002 | -0.010 | 0.005 | 0.747 |
| V4 | class | bray_curtis | M3_rf - M1_knn5 | 24 | 0.004 | -0.005 | 0.013 | 0.565 |
| V4 | class | f1 | M3_rf - M1_knn5 | 24 | -0.041 | -0.075 | -0.009 | 0.043 |
| V4 | class | mae | M3_rf - M1_knn5 | 24 | -0.001 | -0.011 | 0.009 | 0.406 |
| V4 | family | bray_curtis | M3_rf - M1_knn5 | 24 | 0.004 | -0.002 | 0.011 | 0.188 |
| V4 | family | f1 | M3_rf - M1_knn5 | 24 | -0.039 | -0.079 | -0.006 | 0.084 |
| V4 | family | mae | M3_rf - M1_knn5 | 24 | -0.005 | -0.013 | 0.004 | 0.151 |
| V4 | genus | bray_curtis | M3_rf - M1_knn5 | 24 | 0.004 | -0.003 | 0.011 | 0.663 |
| V4 | genus | f1 | M3_rf - M1_knn5 | 24 | -0.052 | -0.091 | -0.020 | 0.040 |
| V4 | genus | mae | M3_rf - M1_knn5 | 24 | -0.007 | -0.016 | 0.002 | 0.151 |
| V4 | order | bray_curtis | M3_rf - M1_knn5 | 24 | 0.004 | -0.002 | 0.011 | 0.264 |
| V4 | order | f1 | M3_rf - M1_knn5 | 24 | -0.035 | -0.067 | -0.007 | 0.101 |
| V4 | order | mae | M3_rf - M1_knn5 | 24 | -0.006 | -0.013 | 0.003 | 0.074 |
| V4 | random | bray_curtis | M3_rf - M1_knn5 | 24 | 0.004 | -0.002 | 0.010 | 0.208 |
| V4 | random | f1 | M3_rf - M1_knn5 | 24 | -0.039 | -0.071 | -0.010 | 0.084 |
| V4 | random | mae | M3_rf - M1_knn5 | 24 | -0.003 | -0.012 | 0.006 | 0.491 |
| full | class | bray_curtis | M3_rf - M1_knn5 | 24 | 0.003 | -0.004 | 0.009 | 0.491 |
| full | class | f1 | M3_rf - M1_knn5 | 24 | -0.049 | -0.083 | -0.019 | 0.004 |
| full | class | mae | M3_rf - M1_knn5 | 24 | -0.005 | -0.012 | 0.001 | 0.095 |
| full | family | bray_curtis | M3_rf - M1_knn5 | 24 | 0.005 | -0.001 | 0.011 | 0.060 |
| full | family | f1 | M3_rf - M1_knn5 | 24 | -0.020 | -0.048 | 0.002 | 0.128 |
| full | family | mae | M3_rf - M1_knn5 | 24 | -0.005 | -0.012 | 0.001 | 0.178 |
| full | genus | bray_curtis | M3_rf - M1_knn5 | 24 | 0.003 | -0.003 | 0.009 | 0.360 |
| full | genus | f1 | M3_rf - M1_knn5 | 24 | -0.052 | -0.084 | -0.024 | 0.002 |
| full | genus | mae | M3_rf - M1_knn5 | 24 | -0.005 | -0.012 | 0.002 | 0.360 |
| full | order | bray_curtis | M3_rf - M1_knn5 | 24 | 0.008 | 0.003 | 0.013 | 0.005 |
| full | order | f1 | M3_rf - M1_knn5 | 24 | -0.036 | -0.068 | -0.012 | 0.015 |
| full | order | mae | M3_rf - M1_knn5 | 24 | -0.004 | -0.010 | 0.003 | 0.240 |
| full | random | bray_curtis | M3_rf - M1_knn5 | 24 | 0.005 | -0.001 | 0.010 | 0.023 |
| full | random | f1 | M3_rf - M1_knn5 | 24 | -0.038 | -0.069 | -0.013 | 0.025 |
| full | random | mae | M3_rf - M1_knn5 | 24 | -0.004 | -0.012 | 0.004 | 0.375 |
| V3V4 | class | bray_curtis | M4_lgbm - M1_knn5 | 24 | 0.002 | -0.009 | 0.013 | 0.966 |
| V3V4 | class | mae | M4_lgbm - M1_knn5 | 24 | -0.009 | -0.022 | 0.004 | 0.019 |
| V3V4 | family | bray_curtis | M4_lgbm - M1_knn5 | 24 | 0.001 | -0.009 | 0.011 | 0.705 |
| V3V4 | family | mae | M4_lgbm - M1_knn5 | 24 | -0.015 | -0.029 | -0.001 | 0.046 |
| V3V4 | genus | bray_curtis | M4_lgbm - M1_knn5 | 24 | 0.002 | -0.008 | 0.012 | 0.546 |
| V3V4 | genus | mae | M4_lgbm - M1_knn5 | 24 | -0.019 | -0.034 | -0.004 | 0.043 |
| V3V4 | order | bray_curtis | M4_lgbm - M1_knn5 | 24 | 0.005 | -0.005 | 0.014 | 0.169 |
| V3V4 | order | mae | M4_lgbm - M1_knn5 | 24 | -0.013 | -0.027 | 0.001 | 0.084 |
| V3V4 | random | bray_curtis | M4_lgbm - M1_knn5 | 24 | 0.004 | -0.006 | 0.013 | 0.360 |
| V3V4 | random | mae | M4_lgbm - M1_knn5 | 24 | -0.016 | -0.029 | -0.005 | 0.043 |
| V4 | class | bray_curtis | M4_lgbm - M1_knn5 | 24 | 0.001 | -0.011 | 0.012 | 0.944 |
| V4 | class | mae | M4_lgbm - M1_knn5 | 24 | -0.011 | -0.022 | 0.003 | 0.037 |
| V4 | family | bray_curtis | M4_lgbm - M1_knn5 | 24 | -0.000 | -0.009 | 0.008 | 0.966 |
| V4 | family | mae | M4_lgbm - M1_knn5 | 24 | -0.015 | -0.028 | -0.001 | 0.034 |
| V4 | genus | bray_curtis | M4_lgbm - M1_knn5 | 24 | -0.003 | -0.013 | 0.007 | 0.439 |
| V4 | genus | mae | M4_lgbm - M1_knn5 | 24 | -0.020 | -0.033 | -0.006 | 0.021 |
| V4 | order | bray_curtis | M4_lgbm - M1_knn5 | 24 | 0.002 | -0.006 | 0.010 | 0.603 |
| V4 | order | mae | M4_lgbm - M1_knn5 | 24 | -0.016 | -0.028 | -0.003 | 0.016 |
| V4 | random | bray_curtis | M4_lgbm - M1_knn5 | 24 | 0.002 | -0.006 | 0.010 | 0.684 |
| V4 | random | mae | M4_lgbm - M1_knn5 | 24 | -0.014 | -0.026 | -0.001 | 0.065 |
| full | class | bray_curtis | M4_lgbm - M1_knn5 | 24 | -0.009 | -0.018 | 0.001 | 0.065 |
| full | class | mae | M4_lgbm - M1_knn5 | 24 | -0.018 | -0.029 | -0.007 | 0.003 |
| full | family | bray_curtis | M4_lgbm - M1_knn5 | 24 | -0.002 | -0.010 | 0.007 | 0.684 |
| full | family | mae | M4_lgbm - M1_knn5 | 24 | -0.019 | -0.033 | -0.006 | 0.016 |
| full | genus | bray_curtis | M4_lgbm - M1_knn5 | 24 | -0.007 | -0.015 | 0.001 | 0.136 |
| full | genus | mae | M4_lgbm - M1_knn5 | 24 | -0.024 | -0.037 | -0.012 | 0.003 |
| full | order | bray_curtis | M4_lgbm - M1_knn5 | 24 | 0.005 | -0.002 | 0.013 | 0.188 |
| full | order | mae | M4_lgbm - M1_knn5 | 24 | -0.016 | -0.027 | -0.004 | 0.014 |
| full | random | bray_curtis | M4_lgbm - M1_knn5 | 24 | 0.000 | -0.007 | 0.008 | 0.747 |
| full | random | mae | M4_lgbm - M1_knn5 | 24 | -0.019 | -0.032 | -0.007 | 0.016 |

## Amplicon arms (random forest gain over null)

![fig3](fig3_amplicon_arms.png)

> **Figure 3. Full-length 16S against in-silico amplicons.** The random forest is trained and tested on the full gene (*full*) or only on the region a primer pair would amplify (in-silico PCR with the primers in `primers.tsv`). Genomes whose 16S does not amplify are left out of that arm. Points are mean gain over the null, per holdout level.
>
> Panels: **KO set F1** = overlap between the predicted and the true set of KEGG orthologs of a genome (1 = perfect); **module completeness MAE** = mean absolute error of predicted KEGG module completeness (0-1 scale; 0 = perfect); **COG profile Bray-Curtis** = dissimilarity between predicted and true COG-category composition (0 = identical).
>
> **How to read it:** the gap between the *full* line and an amplicon line is the information lost by sequencing only that region.
>
> **Here (KO F1 gain):** random: full 0.037, V3V4 0.027, V4 0.025; class: full 0.015, V3V4 0.007, V4 0.007.

| target_type | metric | arm | level | mean_gain |
|---|---|---|---|---|
| ko | f1 | full | random | 0.037 |
| ko | f1 | full | genus | 0.042 |
| ko | f1 | full | family | 0.043 |
| ko | f1 | full | order | 0.025 |
| ko | f1 | full | class | 0.015 |
| ko | f1 | V3V4 | random | 0.027 |
| ko | f1 | V3V4 | genus | 0.038 |
| ko | f1 | V3V4 | family | 0.031 |
| ko | f1 | V3V4 | order | 0.020 |
| ko | f1 | V3V4 | class | 0.007 |
| ko | f1 | V4 | random | 0.025 |
| ko | f1 | V4 | genus | 0.034 |
| ko | f1 | V4 | family | 0.029 |
| ko | f1 | V4 | order | 0.023 |
| ko | f1 | V4 | class | 0.007 |
| module | mae | full | random | 0.014 |
| module | mae | full | genus | 0.018 |
| module | mae | full | family | 0.014 |
| module | mae | full | order | 0.013 |
| module | mae | full | class | 0.013 |
| module | mae | V3V4 | random | 0.013 |
| module | mae | V3V4 | genus | 0.016 |
| module | mae | V3V4 | family | 0.013 |
| module | mae | V3V4 | order | 0.013 |
| module | mae | V3V4 | class | 0.012 |
| module | mae | V4 | random | 0.009 |
| module | mae | V4 | genus | 0.013 |
| module | mae | V4 | family | 0.011 |
| module | mae | V4 | order | 0.011 |
| module | mae | V4 | class | 0.009 |
| cog | bray_curtis | full | random | 0.004 |
| cog | bray_curtis | full | genus | 0.007 |
| cog | bray_curtis | full | family | 0.006 |
| cog | bray_curtis | full | order | 0.004 |
| cog | bray_curtis | full | class | 0.007 |
| cog | bray_curtis | V3V4 | random | 0.004 |
| cog | bray_curtis | V3V4 | genus | 0.006 |
| cog | bray_curtis | V3V4 | family | 0.004 |
| cog | bray_curtis | V3V4 | order | 0.004 |
| cog | bray_curtis | V3V4 | class | 0.006 |
| cog | bray_curtis | V4 | random | 0.002 |
| cog | bray_curtis | V4 | genus | 0.006 |
| cog | bray_curtis | V4 | family | 0.004 |
| cog | bray_curtis | V4 | order | 0.002 |
| cog | bray_curtis | V4 | class | 0.003 |

## KO prediction: median per-KO AUROC (pooled over folds)

For each KO, how well the predicted score ranks genomes that have it above genomes that lack it (0.5 = chance, 1 = perfect); median over KOs. The null predicts one constant per fold, so its pooled AUROC can fall below 0.5 when held-out clades differ from the training prevalence; it is shown only for reference.

| arm | level | M0_null | M1_knn1 | M1_knn5 | M3_rf |
|---|---|---|---|---|---|
| V3V4 | class | 0.226 | 0.500 | 0.507 | 0.398 |
| V3V4 | family | 0.205 | 0.548 | 0.556 | 0.437 |
| V3V4 | genus | 0.213 | 0.556 | 0.580 | 0.489 |
| V3V4 | order | 0.196 | 0.500 | 0.536 | 0.407 |
| V3V4 | random | 0.213 | 0.535 | 0.559 | 0.473 |
| V4 | class | 0.226 | 0.500 | 0.509 | 0.371 |
| V4 | family | 0.205 | 0.500 | 0.556 | 0.432 |
| V4 | genus | 0.213 | 0.531 | 0.589 | 0.462 |
| V4 | order | 0.196 | 0.477 | 0.555 | 0.391 |
| V4 | random | 0.213 | 0.542 | 0.565 | 0.432 |
| full | class | 0.226 | 0.500 | 0.559 | 0.422 |
| full | family | 0.205 | 0.525 | 0.587 | 0.487 |
| full | genus | 0.213 | 0.542 | 0.621 | 0.505 |
| full | order | 0.196 | 0.500 | 0.548 | 0.464 |
| full | random | 0.213 | 0.544 | 0.587 | 0.477 |

## PICRUSt2 (M2) out of the box: external reference only

Its reference genomes overlap the test genomes, so these scores are optimistic (DESIGN §4.1). Our models shown at the random holdout level, the closest comparison. PICRUSt2 places each 16S in its own reference tree and predicts KO copy numbers by hidden-state prediction; a KO counts as predicted when its copy number is >= 1. Rows split the genomes by whether PICRUSt2's closest reference genome is the genome itself (leakage). NSTI = tree distance to the closest reference genome.

| model | f1 | jaccard | nsti | n_genomes |
|---|---|---|---|---|
| M0_null | 0.502 | 0.345 | nan | 24.000 |
| M1_knn1 | 0.561 | 0.402 | nan | 24.000 |
| M1_knn5 | 0.577 | 0.417 | nan | 24.000 |
| M2_picrust2 (self in reference) | 0.956 | 0.918 | 0.023 | 10.000 |
| M2_picrust2 (self not in reference) | 0.810 | 0.694 | 0.390 | 14.000 |
| M3_rf | 0.539 | 0.383 | nan | 24.000 |

Tables behind every figure are in the `fig*_data.tsv` files next to this report.

"""Markdown report of the benchmark (rendered to HTML by render_html.py). Called by bench_eval.py."""
from pathlib import Path

import numpy as np
import pandas as pd

from bench_figures import METHODS, PRIMARY_LEVELS

LABEL = {**{k: v[0] for k, v in METHODS.items()},
         "s16fun_t05": "s16fun, calls at p >= 0.5 (ablation)", "knn5_ef1": "5 nearest, expected-F1 calls (ablation)",
         "oracle5": "ceiling: 5 functionally nearest genomes (oracle)"}
for _m in ["mp", "emp_prob", "pic", "subtree_average"]:
    LABEL[f"picrust2_{_m}"] = f"PICRUSt2 {_m}" + (" (default)" if _m == "mp" else "")
    LABEL[f"picrust2_{_m}_tuned"] = f"PICRUSt2 {_m}, tuned threshold"


def md_table(df, floatfmt="{:.3f}"):
    cols = list(df.columns)

    def fmt(v):
        if isinstance(v, (float, np.floating)):
            return "" if np.isnan(v) else floatfmt.format(v)
        return str(v)
    return "\n".join(["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)] +
                     ["| " + " | ".join(fmt(v) for v in r) + " |" for r in df.itertuples(index=False)])


def quote(*paras):
    return "\n>\n".join("> " + p for p in paras)


def write_report(out, phase, genomes, summary, prim, sec, novt, best_setting, thresholds, res, ordi,
                 stacker_importance, pg, community=None):
    out = Path(out)
    full = summary[summary.arm == "full"]
    n_genomes = pg[(pg.arm == "full")].accession.nunique()
    wins = prim.win.all() if not prim.empty else False
    L = []
    title = "Benchmark: s16fun against PICRUSt2 and other 16S function predictors"
    L += [f"# {title}", "",
          f"Phase: **{phase}** ({'development set: method design and tuning; not the confirmatory claim' if phase == 'dev' else 'confirm set: frozen methods, scored once; this is the pre-registered claim'}). "
          f"{n_genomes:,} held-out genomes from PICRUSt2's own reference (GTDB r214); see BENCHMARK.md for the pre-registration.", ""]

    # executive summary
    L += ["## Executive summary", ""]
    if not prim.empty:
        verdict = ("**s16fun meets the pre-registered co-primary rule at every holdout level.**" if wins else
                   "**The pre-registered co-primary rule is not met at every level.**")
        L += [f"- {verdict} Rule: per-KO AUROC at least +0.05 above PICRUSt2's best setting (CI above 0) and KO F1 not "
              "worse by more than 0.01 (CI lower bound >= -0.01), at genus, family, order and class.", ""]
        for r in prim.itertuples():
            L += [f"  - **{r.level} held out:** per-KO AUROC {r.auroc_s16fun:.3f} vs {r.auroc_picrust2_best:.3f} "
                  f"({r.auroc_diff:+.3f}, CI {r.auroc_ci_lo:+.3f} to {r.auroc_ci_hi:+.3f}; s16fun higher on "
                  f"{r.kos_s16fun_better:.0%} of {r.kos:,} KOs); KO F1 {r.f1_s16fun:.3f} vs {r.f1_picrust2_best:.3f} "
                  f"({r.f1_diff:+.3f}, CI {r.f1_ci_lo:+.3f} to {r.f1_ci_hi:+.3f})."]
        L += [""]
    if not novt.empty:
        d = novt[novt.bin_by == "nsti_bin"].pivot(index="bin", columns="method", values="ko_f1")
        if {"s16fun", "picrust2_best"} <= set(d.columns):
            d = d.dropna()
            diff = (d.s16fun - d.picrust2_best)
            L += [f"- **When:** by PICRUSt2's own novelty score (NSTI), s16fun minus PICRUSt2 ranges from {diff.min():+.3f} "
                  f"(NSTI {diff.idxmin()}) to {diff.max():+.3f} (NSTI {diff.idxmax()}).", ""]
    ceil = full[full.method.isin(["oracle5", "picrust2_best", "s16fun"])].pivot(index="level", columns="method", values="ko_f1")
    if {"oracle5", "picrust2_best", "s16fun"} <= set(ceil.columns):
        ceil = ceil.reindex(PRIMARY_LEVELS).dropna()
        L += ["- **Headroom.** An oracle that copies the 5 reference genomes whose KO sets are truly most similar (information "
              "no 16S method has) reaches KO F1 " + ", ".join(f"{lv} {r.oracle5:.3f}" for lv, r in ceil.iterrows()) +
              ". That is the ceiling for predicting from relatives; PICRUSt2 sits " +
              ", ".join(f"{r.oracle5 - r.picrust2_best:.3f}" for _, r in ceil.iterrows()) + " below it and s16fun " +
              ", ".join(f"{r.oracle5 - r.s16fun:.3f}" for _, r in ceil.iterrows()) + " (genus to class).", ""]
    if community is not None:
        own = community[(community.copies_from == "own 16S copies")]
        own = own[own.arm == ("V4" if "V4" in set(own.arm) else own.arm.iloc[0])]
        for habitat, h in own.groupby("habitat"):
            h = h.set_index("method").bray_curtis_mean
            p2 = h[h.index.str.startswith("picrust2_")]
            if "s16fun" in h and len(p2):
                L += [f"- **Community KO profiles, {habitat}** ({own.arm.iloc[0]} ASVs from held-out genera): Bray-Curtis to the "
                      f"truth {h['s16fun']:.3f} (s16fun) vs {p2.min():.3f} (PICRUSt2 best, {p2.idxmin().removeprefix('picrust2_')}) "
                      f"and {h.get('picrust2_mp', np.nan):.3f} (PICRUSt2 default): {1 - h['s16fun'] / p2.min():.0%} less error "
                      "than PICRUSt2's best."]
        L += [""]
    if not res.empty:
        def get(p, c):
            v = res[res.process.str.contains(p)][c]
            return v.iloc[0] if len(v) else np.nan
        L += [f"- **Cost per task** (one fold of ~5,500 held-out genomes against ~22,000 reference genomes): s16fun "
              f"{get('BENCH_PREDICT', 'median_task_min'):.1f} min, {get('BENCH_PREDICT', 'peak_gb'):.0f} GB; PICRUSt2 "
              f"{get('PICRUST2_FOLD', 'median_task_min'):.0f} min, {get('PICRUST2_FOLD', 'peak_gb'):.0f} GB "
              "(all four HSP settings; one placement).", ""]
    L += ["- Every method sees the same query sequences, the same reference genomes and labels (PICRUSt2's own), "
          "and the same held-out clade; PICRUSt2 gets a per-fold reference built from its own files, so no test "
          "genome is ever in any method's reference.", ""]

    # methods
    L += ["## Methods compared", "", "![methods](fig0_methods.png)", "",
          quote("**Figure 0. The methods and the fairness rules.** Each row is one method family, left to right from the query "
                "16S to KO calls. All of them use the same reference genomes (PICRUSt2's 27,870 GTDB r214 genomes with "
                "PICRUSt2's own KO copy-number tables), minus the clade held out in that fold.",
                "**s16fun** finds the 25 most similar reference 16S sequences, turns them into evidence for each KO (what the "
                "nearest genomes carry, how common the KO is in the nearest genome's genus, family, order, class and phylum, "
                "how similar the hits are, how conserved the KO is within genera), and a gradient-boosted model learned on "
                "dev genomes turns that evidence into a probability. It then calls the KO set that maximises expected F1. "
                "**PICRUSt2** places the sequence in its 16S tree and predicts copy numbers by hidden-state prediction; all "
                "four of its HSP settings are run, and its best one at each level is the comparator."), "",
          "| method | what it is | trained on |", "|---|---|---|",
          "| s16fun (ours) | vsearch hits -> evidence -> boosted stacker -> expected-F1 calls | stacker: dev genomes' held-out predictions (cross-fitted by family) |",
          "| PICRUSt2 2.6.4 | EPA-ng placement + castor HSP; mp (default), emp_prob, pic, subtree_average; thresholds of continuous settings tuned on dev | per-fold reference built from its own files |",
          "| k-mer neural net | embeRNA design: 1-5-mer frequencies -> 2 x 512 ReLU -> sigmoid per KO | each fold's training genomes, same region as the query |",
          "| 5 nearest genomes | identity-weighted average of the 5 most similar 16S | none |",
          "| nearest genome | copy the KOs of the most similar 16S | none |",
          "| Piphillin / Tax4Fun2-like | nearest genome, only if identity >= 97 % | none |",
          "| MicFunPred / PanFP-like | KOs in >= 50 % of the nearest genome's genus, only if identity >= 97 % | none |",
          "| null | KOs present in >= 50 % of reference genomes | none |", ""]

    # primary
    L += ["## 1. Co-primary endpoints", "", "![primary](fig1_primary.png)", "",
          quote("**Figure 1. The co-primary endpoints.** Left: mean per-genome KO F1 of every method (full-length 16S) when "
                "whole species (random), genera, families, orders or classes are held out of the reference; the black dotted "
                "line is a ceiling, not a method (the KOs of the 5 reference genomes whose KO sets are most similar to the "
                "truth). Middle: (a) per-KO AUROC, the mean over KOs of s16fun minus PICRUSt2's best setting: how much better "
                "s16fun ranks genomes that carry a KO above genomes that do not. Right: (b) per-genome KO F1, s16fun minus "
                "PICRUSt2. Intervals: 95 %, from resampling whole genera (genomes of one genus are not independent).",
                "**How to read it:** (a) is met when the point is right of the +0.05 line and the interval stays above 0; (b) is "
                "met when the interval's left end stays right of the -0.01 line. A level counts only if both are met."), "",
          md_table(prim[["level", "genomes", "genera", "auroc_s16fun", "auroc_picrust2_best", "auroc_picrust2_setting",
                         "auroc_diff", "auroc_ci_lo", "auroc_ci_hi", "auroc_superior", "f1_s16fun", "f1_picrust2_best",
                         "f1_diff", "f1_ci_lo", "f1_ci_hi", "f1_noninferior", "win"]]), ""]
    if not best_setting.empty:
        L += ["PICRUSt2's best setting per level for F1 (chosen after seeing these genomes, generous to PICRUSt2): " +
              ", ".join(f"{r.level}/{r.arm}: {LABEL.get(r.picrust2_best_setting, r.picrust2_best_setting)}"
                        for r in best_setting.itertuples() if r.arm == "full") +
              (". Dev-tuned thresholds for its continuous settings: " + ", ".join(f"{k} {v}" for k, v in thresholds.items())
               if thresholds else "") + ".", ""]
    tab = full.pivot(index="method", columns="level", values="ko_f1").reindex(columns=["random"] + PRIMARY_LEVELS)
    tab = tab.loc[[m for m in LABEL if m in tab.index]]
    tab.index = [LABEL[m] for m in tab.index]
    L += ["All methods, mean KO F1 per genome (full length):", "", md_table(tab.reset_index().rename(columns={"index": "method"})), ""]

    # secondary
    L += ["## 2. Secondary endpoints", "",
          "Paired per genome, s16fun minus PICRUSt2's best setting, signed so that **positive = s16fun better** "
          "(errors sign-flipped); 95 % CI by resampling genera.", ""]
    if not sec.empty:
        s = sec[(sec.arm == "full")].drop(columns=[c for c in ["arm"] if c in sec])
        L += [md_table(s), ""]

    # when and why
    L += ["## 3. When and why each method wins", "", "![novelty](fig2_novelty.png)", "",
          quote("**Figure 2. Accuracy by novelty.** Every held-out genome (all levels pooled, full length) binned by how far it "
                "is from the reference: left, 16S identity to the most similar reference genome; right, PICRUSt2's NSTI. "
                "Points: mean KO F1 per bin.",
                "**How to read it:** where the blue line is above the pink one, s16fun is more accurate for genomes of that "
                "novelty. NSTI is PICRUSt2's own yardstick, so the right panel is the fairest view of where PICRUSt2 is "
                "strong.")]
    L += ["", md_table(novt.pivot_table(index=["bin_by", "bin"], columns="method", values="ko_f1", observed=True)
                       .reset_index()), ""]
    L += ["![phyla](fig3_phyla.png)", "",
          quote("**Figure 3. By phylum** (genus held out, full length, phyla with >= 30 held-out genomes). Grey bars join "
                "the two methods' mean KO F1 for that phylum."), ""]
    L += ["![per KO](fig4_per_ko.png)", "",
          quote("**Figure 4. Ranking genomes per KO** (genus held out). Each dot is a KO: AUROC of the method's score "
                "(s16fun: probability; PICRUSt2: predicted copy number) for separating genomes that carry the KO from "
                "genomes that do not. Above the diagonal: s16fun ranks better."), ""]
    if stacker_importance and Path(stacker_importance).exists():
        imp = pd.read_csv(stacker_importance, sep="\t", index_col=0).gain_share.sort_values(ascending=False)
        L += ["**Why s16fun works (what the stacker relies on)**, share of total split gain: " +
              ", ".join(f"{k} {v:.0%}" for k, v in imp.head(6).items()) + ".", ""]
    abl = full[full.method.isin(["s16fun", "s16fun_t05", "knn5_ef1", "knn5"])].pivot(index="method", columns="level", values="ko_f1")
    if len(abl) > 1:
        abl.index = [LABEL[m] for m in abl.index]
        L += ["Ablations (mean KO F1, full length): stacker vs plain 5-NN, and expected-F1 calls vs a 0.5 cut-off.", "",
              md_table(abl.reindex(columns=["random"] + PRIMARY_LEVELS).reset_index().rename(columns={"index": "variant"})), ""]

    # arms
    L += ["## 4. Amplicons", "", "![arms](fig5_arms.png)", "",
          quote("**Figure 5. Full-length 16S against the V3-V4 and V4 amplicons** (in-silico PCR of the same genomes). "
                "Every method gets the amplicon; references stay full length (PICRUSt2 places the fragment in its full-length "
                "tree; s16fun aligns the fragment to full-length references; the k-mer net is trained on the region)."), ""]

    # ordination
    if ordi is not None:
        st = ordi["stats"].copy()
        st["method"] = st.method.map(lambda m: LABEL.get(m, m))
        L += ["## 5. Predicted against true functional profiles", "", "![ordination](fig6_ordination.png)", "",
              quote("**Figure 6. Do predicted profiles reproduce the functional landscape?** Genus held out, full length, a random "
                    "subset of held-out genomes. Per row, the true profiles (grey) and each method's predictions share one PCoA "
                    "(top: KO presence, Jaccard; bottom: module completeness, Bray-Curtis). Thin lines join a genome's true and "
                    "predicted profile.",
                    "**Numbers:** Mantel r = correlation of true and predicted genome-to-genome distances (1 = which genomes are "
                    "alike is preserved); Procrustes m² = shape mismatch (0 = same); PERMANOVA R² and ANOSIM R = how separable "
                    "the predicted cloud is from the true cloud (0 = not); spread = PERMDISP dispersion relative to truth "
                    "(below 1 = predictions pulled towards typical profiles). 999 permutations swap labels within each genome; "
                    "with ~1,000 paired genomes tiny effects reach p = 0.001, so read the effect sizes."), "",
              md_table(st.drop(columns=["genomes"]), "{:.4f}"), ""]

    # communities
    if community is not None:
        c = community.copy()
        c["method"] = c.method.map(lambda x: LABEL.get(x, x.replace("picrust2_", "PICRUSt2 ")))
        L += ["## 6. Simulated communities", "", "![communities](fig8_community.png)", "",
              quote("**Figure 8. Community KO profiles from ASVs** (the task PICRUSt2 is built for). Communities are drawn only "
                    "from genomes whose genus was held out, one V4 ASV per genome, in three habitat designs (gut, marine, soil: "
                    "ASVs per sample and evenness as in typical surveys). Reads = cells x true 16S copies; each method converts "
                    "reads back to cells with its own 16S copy prediction and to KO abundance with its own copy predictions. "
                    "Boxes: replicate communities (100 samples each).",
                    "**How to read it:** lower is closer to the true KO profile. Solid blue: s16fun as shipped; pale blue: s16fun "
                    "with hard 0/1 calls; pink: PICRUSt2's four settings."), ""]
        for arm in sorted(set(c.arm)):
            for src in ["own 16S copies", "true 16S copies"]:
                t = c[(c.arm == arm) & (c.copies_from == src)].sort_values(["habitat", "bray_curtis_mean"])
                cols = ["habitat", "method", "bray_curtis_mean", "bray_curtis_lo", "bray_curtis_hi", "spearman_mean"] + (
                    ["lfc_spearman_mean", "da_precision_mean", "da_recall_mean"] if src == "own 16S copies" else [])
                L += [f"{arm} ASVs, {src}:", "", md_table(t[cols], "{:.4f}"), ""]

    # resources
    if not res.empty:
        L += ["## 7. Run time and memory", "", "![resources](fig7_resources.png)", "",
              quote("**Figure 7. Cost per benchmark task**, from the Nextflow trace: median wall time and peak memory of one task "
                    "(one holdout fold, all its held-out genomes)."), "", md_table(res), ""]

    # fairness notes
    L += ["## 8. Fairness: what favours whom", "",
          "- **Same data.** Reference genomes, 16S sequences and KO labels are PICRUSt2's own (PICRUSt2-SC, GTDB r214, "
          "eggNOG-mapper annotation). Truth for a held-out genome = its row in PICRUSt2's KO table.",
          "- **No leakage.** Every fold removes the held-out genomes from every reference file (PICRUSt2: alignment, tree, "
          "HMM rebuilt, trait tables; s16fun: sequences and KO matrix); asserted per task. Clade levels hold out whole taxa.",
          "- **In PICRUSt2's favour:** its tree topology and RAxML model parameters were estimated with the held-out genomes "
          "present (pruning removes them as tips, not their influence); the comparator is its best of 8 settings per level, "
          "chosen after seeing the scores.",
          "- **In s16fun's favour:** its stacker was trained on dev genomes of this same reference (cross-fitted by family in "
          "the dev report; trained on all dev genomes and frozen for the confirm report).",
          "- **KOs PICRUSt2 cannot predict** (absent from a fold's reference) count as absent for every method.", ""]

    L += ["## 9. Outputs", "",
          "- `per_genome.parquet`: every genome x method x level x arm: KO F1, precision, recall, Jaccard, module and copy-number error.",
          "- `summary.tsv`, `primary.tsv`, `secondary.tsv`, `by_novelty.tsv`, `per_ko_auroc.tsv.gz`, `ordination_stats.tsv`, "
          "`resources.tsv`, `picrust2_best_setting.tsv`, `picrust2_thresholds.tsv` (dev).", ""]
    (out / "benchmark_report.md").write_text("\n".join(L) + "\n")

#!/usr/bin/env python3
"""Exploration report: every idea tried, its rationale, its result against PICRUSt2, and the figures.

Usage: python exploration/report.py --results results/exploration/results.tsv --outdir reports/exploration
Writes exploration_report.md (+ figures); render to HTML with bin/render_html.py.
"""
import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from bench_figures import GRID, INK, INK2, style  # noqa: E402
from bench_report import md_table, quote  # noqa: E402

LEVELS = ["random", "genus", "family", "order", "class"]
BLUE, PINK, GREY = "#2a78d6", "#d55181", "#8a8984"

IDEAS = [  # id, idea, rationale
    ("E1", "5-NN copy averaging", "Baseline: copy the KO copy numbers and 16S copies of the 5 most similar references."),
    ("E2", "s16fun: probability x copies when present", "Expected value instead of a single best guess: the community "
     "profile is a sum over ASVs, and the mean minimises squared error while a mode is biased for uncertain KOs."),
    ("E3", "E2 + learned 16S copy model", "16S copy normalisation is half of the community error (D1); a model on the copy "
     "numbers and identities of the nearest references and clade means."),
    ("E3h", "E3 with hard 0/1 calls", "Ablation: isolates the value of expectations over calls."),
    ("E4", "direct copy regression (Tweedie boosting)", "Optimise E[copies] directly instead of p x E[copies | present]."),
    ("E5", "ratio model: KO copies per 16S copy", "A community profile only needs, per ASV, KO copies divided by 16S copies. "
     "Estimating the ratio from relatives avoids dividing two noisy estimates whose errors are correlated."),
    ("E5b", "ratio model with the stacker's probability", "E5 restricted to relatives carrying the KO, times the presence probability."),
    ("E6", "learned ratio: KO copies per 16S copy (Tweedie boosting)", "E4 and E5 combined: regress the ratio a read "
     "carries directly on the ratio and presence evidence of relatives, so no 16S copy estimate is ever divided."),
    ("D1", "diagnostic: true 16S copies", "Upper bound if 16S copy number were known."),
    ("D2", "diagnostic: s16fun KOs + PICRUSt2 16S copies", "Is PICRUSt2's 16S copy prediction the missing piece?"),
    ("D3", "diagnostic: mean of s16fun and PICRUSt2 pic", "Do the two methods make complementary errors?"),
]
PROTOTYPE = [  # earlier dev-fold prototypes (genus fold 0, full length), D086/D093
    ("Low-rank (SVD) denoising of the 5/25-NN profile as stacker features", "KO F1 +0.001", "dropped"),
    ("Consensus taxon of the top 10 hits (instead of the top hit) for clade priors", "KO F1 +0.001", "dropped"),
    ("Phylogenetic kernel on PICRUSt2's reference tree (patristic distance from the nearest hit)", "KO F1 -0.0005", "dropped"),
    ("Tree-kernel and learned estimators for 16S copy number", "log error 0.46-0.50 vs PICRUSt2 0.43", "learned model kept (best of ours)"),
    ("Ensemble of s16fun and PICRUSt2 calls", "KO F1 +0.005 (0.836 vs 0.831)", "not a standalone method"),
    ("Uncertainty-filtered differential abundance", "precision 0.046 vs 0.085", "dropped"),
    ("Expected-F1 calls instead of a 0.5 cut-off", "KO F1 +0.004", "kept"),
]


def short(v):
    return v.split(" ")[0] if v[:1] in "ED" and v[1:2].isdigit() else v


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    r = pd.read_csv(a.results, sep="\t")
    com = r[r.metric == "community_bray_curtis"]
    variants = list(dict.fromkeys(r.variant))
    p2 = [v for v in variants if v.startswith("PICRUSt2")]

    # community BC per variant, arm, habitat (mean over levels and replicates)
    tab = com.groupby(["arm", "habitat", "variant"]).value.mean().unstack("variant")
    best_p2 = tab[p2].min(axis=1)
    rel = tab.div(best_p2, axis=0) - 1          # relative to PICRUSt2's best setting (negative = better)
    rel.to_csv(out / "community_vs_picrust2.tsv", sep="\t")

    # figure 1: community error per variant, V4 ASVs, by habitat
    arm = "V4" if "V4" in tab.index.get_level_values(0) else tab.index.get_level_values(0)[0]
    t = tab.loc[arm]
    order = t.mean().sort_values().index.tolist()
    fig, ax = plt.subplots(figsize=(11, 0.42 * len(order) + 1.5))
    y = np.arange(len(order))
    for j, (h, mk) in enumerate(zip(t.index, ["o", "s", "D"])):
        ax.scatter(t.loc[h, order], y + (j - 1) * 0.18, marker=mk, s=36,
                   color=[PINK if v.startswith("PICRUSt2") else GREY if v.startswith("D") else BLUE for v in order],
                   label=h, zorder=3)
    ax.axvline(best_p2.loc[arm].min(), color=PINK, lw=1, ls=":")
    ax.set_yticks(y, order, fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel(f"Bray-Curtis to the true KO profile ({arm} ASVs; lower is better)", color=INK2)
    style(ax)
    ax.grid(axis="x", color=GRID)
    ax.legend(title="habitat (marker)", frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "fig1_community_variants.png", dpi=150)
    plt.close(fig)

    # figure 2: by holdout level, main variants
    keyv = [v for v in ["E3 s16fun p x copies, learned 16S", "E5b ratio, p x ratio when present",
                        "E4 copy regression (Tweedie), learned 16S", "E6 ratio regression (Tweedie)"] if v in variants]
    lv = com[com.arm == arm].groupby(["level", "variant"]).value.mean().unstack()
    fig, ax = plt.subplots(figsize=(8, 4.2))
    x = np.arange(len(LEVELS))
    ax.plot(x, lv.reindex(LEVELS)[p2].min(axis=1), marker="s", color=PINK, label="PICRUSt2 (best setting per level)")
    for v, mk in zip(keyv, ["o", "^", "D", "v"]):
        ax.plot(x, lv.reindex(LEVELS)[v], marker=mk, color=BLUE, alpha=1 if v.startswith("E3") else 0.6, label=v)
    ax.set_xticks(x, ["random split" if lv_ == "random" else f"{lv_} held out" for lv_ in LEVELS])
    ax.set_ylabel(f"Bray-Curtis to truth ({arm}, all habitats)", color=INK2)
    style(ax)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "fig2_by_level.png", dpi=150)
    plt.close(fig)

    # per-genome metrics table
    pg = r[r.metric.isin(["ko_f1", "ko_auroc_mean", "copy_mae", "ssu_log_error"])]
    pgt = pg[pg.arm == arm].groupby(["variant", "metric"]).value.mean().unstack().reindex(variants).dropna(how="all")
    lvl = com[com.arm == arm].groupby(["variant", "level"]).value.mean().unstack().reindex(columns=LEVELS)
    da = r[(r.metric == "community_lfc_spearman") & (r.arm == arm)].groupby("variant").value.mean()

    s16 = "E3 s16fun p x copies, learned 16S"
    L = ["# Exploration: beating PICRUSt2 on community function profiles from ASVs", "",
         "Every idea tried before the benchmark freeze, with its rationale and its result against PICRUSt2. Dev genomes only "
         "(the confirm genomes are untouched); references exclude the held-out clade; every learned piece is cross-fitted by "
         "family. Communities: built only from held-out genomes, one ASV per genome, gut / marine / soil designs.", "",
         "## Executive summary", ""]
    if s16 in tab:
        for (arm_, h), v in rel[s16].items():
            if arm_ == arm:
                L += [f"- **{h}:** s16fun (E3) community error {tab.loc[(arm_, h), s16]:.3f} vs PICRUSt2's best "
                      f"{best_p2.loc[(arm_, h)]:.3f}: **{-v:.0%} {'less' if v < 0 else 'more'}**."]
    winners = rel.loc[arm].mean().sort_values()
    L += [f"- Best variant overall ({arm} ASVs, mean over habitats): **{winners.index[0]}** "
          f"({-winners.iloc[0]:.0%} less error than PICRUSt2's best setting).", ""]

    L += ["## Ideas tested", "", "| id | idea | why it could work |", "|---|---|---|"]
    L += [f"| {i} | {n} | {why} |" for i, n, why in IDEAS]
    L += ["", "## Result: community profiles", "", "![variants](fig1_community_variants.png)", "",
          quote(f"**Figure 1. Community KO profile error of every variant** ({arm} ASVs, mean over holdout levels and "
                "replicate communities). Each marker is one habitat design; pink = PICRUSt2's settings, blue = our "
                "variants, grey = diagnostics that use information a user would not have (true 16S copies) or PICRUSt2's "
                "output. The dotted line is PICRUSt2's best setting.",
                "**How to read it:** further left is closer to the true community KO profile."), "",
          "Error relative to PICRUSt2's best setting (negative = less error than PICRUSt2):", "",
          md_table((rel * 100).round(1).reset_index(), "{:.1f}"), "",
          "![levels](fig2_by_level.png)", "",
          quote("**Figure 2. By holdout level.** Community error when the ASVs' whole genus, family, order or class is "
                "missing from the reference."), "",
          "Community error by holdout level:", "", md_table(lvl.reset_index(), "{:.4f}"), "",
          "Differential-abundance agreement (Spearman of predicted vs true log2 fold changes, mean):", "",
          md_table(da.rename("lfc_spearman").reset_index(), "{:.3f}"), "",
          "## Result: per genome", "", md_table(pgt.reset_index(), "{:.4f}"), "",
          "## Earlier prototypes (genus fold 0, full length)", "", "| idea | result | verdict |", "|---|---|---|"]
    L += [f"| {a_} | {b} | {c} |" for a_, b, c in PROTOTYPE]
    L += ["", "## Caveats", "",
          "- Dev folds only; the numbers that count are the frozen confirm run.",
          "- Several ideas were tried on the same dev genomes (selection); the confirm lockbox exists for this reason.",
          "- One ASV per genome; intragenomic 16S variants and sequencing error are not simulated.",
          "- Habitat designs are parametric (DESIGN §B.5)."]
    (out / "exploration_report.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()

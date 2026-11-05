#!/usr/bin/env python3
"""Aggregate Aim 1 results (all arm x level task folders) into tables, figures and a markdown report.

Usage: aim1_report.py --runs aim1_dirs... --outdir aim1_report [--seed 42] [--boot 1000]

Headline per-genome metrics: KO F1 (higher better), module MAE (lower), COG Bray-Curtis (lower).
"Gain" = improvement over the null model on the same genome, signed so that positive = better.
Paired comparisons (learned model vs k-NN) use per-genome differences: mean, bootstrap 95 % CI
over genomes, and the Wilcoxon signed-rank p-value.
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import wilcoxon  # noqa: E402

HIGHER_BETTER = {"f1": True, "jaccard": True, "spearman": True, "mae": False, "bray_curtis": False}
HEADLINE = [("ko", "f1", "KO set F1"), ("module", "mae", "Module completeness MAE"),
            ("cog", "bray_curtis", "COG profile Bray-Curtis")]
LEVELS = ["random", "genus", "family", "order", "class"]
# Fixed categorical order (validated: dataviz palette slots 1-4); null is a neutral dashed baseline.
MODEL_STYLE = {
    "M0_null": dict(color="#8a8984", ls="--", label="null"),
    "M1_knn1": dict(color="#2a78d6", ls="-", label="1-NN"),
    "M1_knn5": dict(color="#eb6834", ls="-", label="5-NN"),
    "M3_rf": dict(color="#1baf7a", ls="-", label="random forest"),
    "M4_lgbm": dict(color="#eda100", ls="-", label="LightGBM"),
}
ARM_COLOR = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
ID_BINS = [0, 0.85, 0.90, 0.94, 0.97, 0.99, 1.0001]
ID_LABELS = ["<85", "85-90", "90-94", "94-97", "97-99", ">=99"]
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def style():
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
        "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False, "font.size": 10, "lines.linewidth": 2,
        "lines.markersize": 6, "legend.frameon": False,
    })


def boot_ci(x, rng, n):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    if len(x) < 2:
        return np.nan, np.nan
    means = rng.choice(x, size=(n, len(x)), replace=True).mean(1)
    return np.percentile(means, 2.5), np.percentile(means, 97.5)


def signed(df):
    """Value with sign flipped for lower-is-better metrics, so larger is always better."""
    return np.where(df["metric"].map(HIGHER_BETTER), df["value"], -df["value"])


def gains(g):
    """Per-genome improvement over the null model (positive = better)."""
    key = ["arm", "level", "accession", "target_type", "metric"]
    s = g.assign(s=signed(g))
    null = s[s.model == "M0_null"].set_index(key)["s"]
    out = s[s.model != "M0_null"].join(null.rename("null_s"), on=key)
    return out.assign(gain=out["s"] - out["null_s"])


def paired(g, a, b, rng, n_boot):
    key = ["arm", "level", "accession", "target_type", "metric"]
    s = g.assign(s=signed(g)).pivot_table(index=key, columns="model", values="s").dropna(subset=[a, b])
    rows = []
    for (arm, level, tt, metric), d in s.groupby(level=["arm", "level", "target_type", "metric"]):
        diff = (d[a] - d[b]).to_numpy()
        lo, hi = boot_ci(diff, rng, n_boot)
        p = wilcoxon(diff).pvalue if np.any(diff != 0) and len(diff) > 5 else np.nan
        rows.append({"arm": arm, "level": level, "target_type": tt, "metric": metric, "comparison": f"{a} - {b}",
                     "n": len(diff), "mean_diff": diff.mean(), "ci_lo": lo, "ci_hi": hi, "wilcoxon_p": p})
    return pd.DataFrame(rows)


def fig_legend(fig, axes):
    """One legend for all panels (a panel may lack a model, e.g. LightGBM has no KO panel)."""
    seen = {}
    for ax in axes:
        for h, l in zip(*ax.get_legend_handles_labels()):
            seen.setdefault(l, h)
    fig.legend(seen.values(), seen.keys(), loc="outside upper center", ncol=len(seen), fontsize=9)


def fig_distance(g, near, out, col="nn_identity", bins=ID_BINS, labels=ID_LABELS,
                 xlabel="16S identity to nearest training genome (%)", fname="fig1_accuracy_vs_identity.png"):
    """Mean per-genome metric by distance-to-training bin, one line per model (full-length arm)."""
    d = g[g.arm == "full"].merge(near[near.arm == "full"], on=["arm", "level", "accession", "fold"])
    d["bin"] = pd.cut(d[col], bins, labels=labels, right=False)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    table = []
    for ax, (tt, metric, title) in zip(axes, HEADLINE):
        sub = d[(d.target_type == tt) & (d.metric == metric)]
        m = sub.groupby(["model", "bin"], observed=False)["value"].agg(["mean", "size"]).reset_index()
        for model, st in MODEL_STYLE.items():
            mm = m[m.model == model]
            if not mm.empty:
                ax.plot(np.arange(len(labels)), mm["mean"], ls=st["ls"], color=st["color"], marker="o",
                        label=st["label"])
        counts = m[m.model == "M0_null"]["size"].tolist()
        ax.set_xticks(range(len(labels)), [f"{l}\nn={c}" for l, c in zip(labels, counts)], fontsize=8)
        ax.set_title(f"{title} ({'higher' if HIGHER_BETTER[metric] else 'lower'} is better)", fontsize=10)
        ax.set_xlabel(xlabel)
        table.append(m.assign(target_type=tt, metric=metric))
    axes[0].set_ylabel("mean per-genome value")
    fig_legend(fig, axes)
    fig.savefig(out / fname, dpi=150)
    plt.close(fig)
    return pd.concat(table)


def nearest_patristic(D, splits):
    """Per level and test genome: patristic distance to the closest genome in the other folds."""
    rows = []
    for level, s in splits.groupby("level"):
        s = s[s.accession.isin(D.index)]
        for f, te in s.groupby("fold"):
            tr = s.loc[s.fold != f, "accession"]
            sub = D.loc[te.accession, tr].to_numpy()
            with np.errstate(all="ignore"):
                nn = np.where(np.isnan(sub).all(1), np.nan, np.nanmin(np.where(np.isnan(sub), np.inf, sub), 1))
            rows.append(pd.DataFrame({"level": level, "fold": f, "accession": te.accession.to_numpy(),
                                      "nn_patristic": nn}))
    return pd.concat(rows)


def fig_gain_by_level(gn, rng, n_boot, out):
    d = gn[gn.arm == "full"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    rows = []
    levels = [l for l in LEVELS if l in d.level.unique()]
    for ax, (tt, metric, title) in zip(axes, HEADLINE):
        models = [m for m in MODEL_STYLE if m != "M0_null" and m in d[d.target_type == tt].model.unique()]
        for i, model in enumerate(models):
            st = MODEL_STYLE[model]
            ys, lo, hi = [], [], []
            for level in levels:
                v = d[(d.target_type == tt) & (d.metric == metric) & (d.model == model) & (d.level == level)]["gain"]
                c = boot_ci(v, rng, n_boot)
                ys.append(v.mean()); lo.append(c[0]); hi.append(c[1])
                rows.append({"target_type": tt, "metric": metric, "model": model, "level": level,
                             "mean_gain": v.mean(), "ci_lo": c[0], "ci_hi": c[1], "n": v.notna().sum()})
            x = np.arange(len(levels)) + (i - (len(models) - 1) / 2) * 0.15
            ax.errorbar(x, ys, yerr=[np.subtract(ys, lo), np.subtract(hi, ys)], fmt="o", color=st["color"],
                        label=st["label"], capsize=0, elinewidth=1.5)
        ax.axhline(0, color=INK2, lw=1)
        ax.set_xticks(range(len(levels)), levels)
        ax.set_xlabel("holdout level (taxa never seen in training)")
        ax.set_title(f"{title}: gain over null", fontsize=10)
    axes[0].set_ylabel("mean per-genome gain (positive = better)")
    fig_legend(fig, axes)
    fig.savefig(out / "fig2_gain_over_null_by_level.png", dpi=150)
    plt.close(fig)
    return pd.DataFrame(rows)


def fig_arms(gn, out, model="M3_rf"):
    arms = ["full"] + sorted(a for a in gn.arm.unique() if a != "full")
    levels = [l for l in LEVELS if l in gn.level.unique()]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    rows = []
    for ax, (tt, metric, title) in zip(axes, HEADLINE):
        for i, arm in enumerate(arms):
            d = gn[(gn.arm == arm) & (gn.model == model) & (gn.target_type == tt) & (gn.metric == metric)]
            m = d.groupby("level")["gain"].mean().reindex(levels)
            ax.plot(range(len(levels)), m.values, marker="o", color=ARM_COLOR[i % len(ARM_COLOR)], label=arm)
            rows += [{"target_type": tt, "metric": metric, "arm": arm, "level": l, "mean_gain": v} for l, v in m.items()]
        ax.axhline(0, color=INK2, lw=1)
        ax.set_xticks(range(len(levels)), levels)
        ax.set_title(f"{title}: {MODEL_STYLE[model]['label']} gain over null", fontsize=10)
        ax.set_xlabel("holdout level")
    axes[0].set_ylabel("mean per-genome gain (positive = better)")
    fig_legend(fig, axes)
    fig.savefig(out / "fig3_amplicon_arms.png", dpi=150)
    plt.close(fig)
    return pd.DataFrame(rows)


def md_table(df, floatfmt=".3f"):
    cols = list(df.columns)
    fmt = lambda v: f"{v:{floatfmt}}" if isinstance(v, (float, np.floating)) else str(v)  # noqa: E731
    return "\n".join(["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
                     + ["| " + " | ".join(fmt(v) for v in r) + " |" for r in df.itertuples(index=False)])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--patristic", help="patristic.parquet (optional)")
    ap.add_argument("--splits", help="splits.tsv (needed with --patristic)")
    ap.add_argument("--picrust2", help="PICRUSt2 genome metrics (optional; external reference, D-log M2)")
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    style()

    g = pd.concat(pd.read_csv(Path(r) / "genome_metrics.tsv", sep="\t") for r in a.runs)
    t = pd.concat(pd.read_csv(Path(r) / "target_metrics.tsv", sep="\t") for r in a.runs)
    near = pd.concat(pd.read_csv(Path(r) / "nearest.tsv", sep="\t") for r in a.runs)

    summary = (g.groupby(["arm", "level", "target_type", "metric", "model"])["value"]
               .agg(["mean", "median", "size"]).reset_index())
    summary.to_csv(out / "summary_genome_metrics.tsv", sep="\t", index=False)
    gn = gains(g)
    comps = pd.concat([paired(g, "M3_rf", "M1_knn1", rng, a.boot), paired(g, "M3_rf", "M1_knn5", rng, a.boot),
                       paired(g[g.target_type != "ko"], "M4_lgbm", "M1_knn5", rng, a.boot)])
    comps.to_csv(out / "paired_comparisons.tsv", sep="\t", index=False)
    curve = fig_distance(g, near, out)
    curve.to_csv(out / "fig1_data.tsv", sep="\t", index=False)
    phylo_lines = []
    if a.patristic:
        npd = nearest_patristic(pd.read_parquet(a.patristic), pd.read_csv(a.splits, sep="\t"))
        near = near.merge(npd, on=["level", "fold", "accession"], how="left")
        qs = np.nanquantile(near.loc[near.arm == "full", "nn_patristic"], [0, .2, .4, .6, .8, 1])
        qs[-1] += 1e-6
        labels = [f"{lo:.2f}-{hi:.2f}" for lo, hi in zip(qs[:-1], qs[1:])]
        fig_distance(g, near, out, "nn_patristic", qs, labels, "GTDB-tree distance to nearest training genome "
                     "(quintiles)", "fig4_accuracy_vs_tree_distance.png").to_csv(out / "fig4_data.tsv", sep="\t", index=False)
        full = near[near.arm == "full"].dropna(subset=["nn_patristic"])
        rho = full[["nn_identity", "nn_patristic"]].corr(method="spearman").iloc[0, 1]
        phylo_lines = ["## Accuracy vs GTDB-tree distance (full length)", "",
                       f"Spearman correlation between 16S identity and tree distance to the nearest training "
                       f"genome: {rho:.3f}.", "", "![fig4](fig4_accuracy_vs_tree_distance.png)", ""]
    pc_lines = []
    if a.picrust2:
        pc = pd.read_csv(a.picrust2, sep="\t")
        self_ref = pc[pc.metric == "self_in_reference"].set_index("accession")["value"] == 1
        pc = pc[pc.metric.isin(["f1", "jaccard", "nsti"])].copy()
        pc["model"] = np.where(pc.accession.map(self_ref), "M2_picrust2 (self in reference)",
                               "M2_picrust2 (self not in reference)")
        ours = g[(g.arm == "full") & (g.level == "random") & (g.target_type == "ko") & g.metric.isin(["f1", "jaccard"])]
        both = pd.concat([pc, ours]).groupby(["model", "metric"])["value"].agg(["mean", "size"])
        both = both["mean"].unstack("metric").join(both["size"].unstack("metric").iloc[:, 0].rename("n_genomes"))
        both = both.reset_index()
        pc_lines = ["## PICRUSt2 (M2) out of the box: external reference only", "",
                    "Its reference genomes overlap the test genomes, so these scores are optimistic (DESIGN §4.1). "
                    "Our models shown at the random holdout level, the closest comparison.", "",
                    md_table(both.round(3)), ""]
    by_level = fig_gain_by_level(gn, rng, a.boot, out)
    by_level.to_csv(out / "fig2_data.tsv", sep="\t", index=False)
    arms = fig_arms(gn, out)
    arms.to_csv(out / "fig3_data.tsv", sep="\t", index=False)
    auroc = (t[(t.target_type == "ko") & (t.metric == "auroc")].groupby(["arm", "level", "model"])["value"]
             .median().unstack("model").reset_index())
    auroc.to_csv(out / "ko_auroc_median.tsv", sep="\t", index=False)
    nn = near[near.arm == "full"].groupby("level")["nn_identity"].describe(percentiles=[.1, .5, .9]).reset_index()

    head = [(tt, m) for tt, m, _ in HEADLINE]
    lines = [
        "# Aim 1 report: predicting function from 16S", "",
        f"Genomes: {g.accession.nunique()}. Arms: {', '.join(sorted(g.arm.unique()))}. "
        f"Levels: {', '.join(l for l in LEVELS if l in g.level.unique())}.", "",
        "## 16S identity to the nearest training genome (full length)", "",
        md_table(nn[["level", "count", "10%", "50%", "90%"]]), "",
        "## Accuracy vs 16S identity (full length, all levels pooled)", "",
        "![fig1](fig1_accuracy_vs_identity.png)", "",
        *phylo_lines,
        "## Gain over the null model by holdout level (full length)", "",
        "![fig2](fig2_gain_over_null_by_level.png)", "",
        md_table(by_level.round(4)), "",
        "## Learned models vs k-NN (paired per genome; positive = learned model better)", "",
        md_table(comps[comps.apply(lambda r: (r.target_type, r.metric) in head, axis=1)]
                 .drop(columns=["target_type"]).round(4)), "",
        "## Amplicon arms (random forest gain over null)", "",
        "![fig3](fig3_amplicon_arms.png)", "",
        md_table(arms.round(4)), "",
        "## KO prediction: median per-KO AUROC (pooled over folds)", "",
        md_table(auroc.round(3)), "",
        *pc_lines,
        "Tables behind every figure are in the `fig*_data.tsv` files next to this report.",
    ]
    (out / "aim1_report.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

"""Benchmark figures (matplotlib PNG). Colours: validated categorical palette (dataviz), fixed per method."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

INK, INK2, GRID, TRUTH = "#1f1e1c", "#52514e", "#e4e3df", "#b9b8b3"
# key: (label, colour, marker). Order = legend order. Validated: adjacent CVD dE >= 9.2 (light).
METHODS = {
    "s16fun": ("s16fun (ours)", "#2a78d6", "o"),
    "picrust2_best": ("PICRUSt2 (best setting)", "#d55181", "s"),
    "kmer_mlp": ("k-mer neural net (embeRNA-like)", "#008300", "^"),
    "knn5": ("5 nearest genomes", "#eda100", "D"),
    "nn1": ("nearest genome", "#6250d6", "v"),
    "nearest97": ("nearest genome if >= 97 % (Piphillin / Tax4Fun2-like)", "#eb6834", "P"),
    "genuscore97": ("genus core if >= 97 % (MicFunPred / PanFP-like)", "#1baf7a", "X"),
    "null": ("null (KO prevalence)", "#8a8984", "_"),
}
PRIMARY_LEVELS = ["genus", "family", "order", "class"]
LEVELS = ["random", "genus", "family", "order", "class"]


def style(ax):
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9)


def fig_methods_diagram(path):
    """What each method does with a 16S sequence, and how the benchmark keeps the comparison fair."""
    fig, ax = plt.subplots(figsize=(16, 7.8))
    ax.set_xlim(0, 16.2)
    ax.set_ylim(0, 7.8)
    ax.axis("off")

    def box(x, y, w, h, text, color, fs=8.5, bold=False):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.03,rounding_size=0.1",
                                    fc=color + "22", ec=color, lw=1.3))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=INK,
                fontweight="bold" if bold else "normal")

    def arrow(x0, y0, x1, y1, color=INK2, ls="-"):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=11, color=color, lw=1.1,
                                     linestyle=ls))

    ax.text(8.1, 7.5, "Same query sequences, same reference genomes and labels, same held-out clade for every method",
            ha="center", fontsize=11.5, color=INK, fontweight="bold")
    box(0.15, 5.2, 2.3, 1.5, "query 16S\n(held-out genome;\nfull length, V4\nor V3-V4)", INK2, 9, True)
    box(0.15, 0.8, 2.3, 2.6, "reference\n= training genomes\n\nPICRUSt2's own 27,870\nGTDB r214 genomes\nminus the held-out\nclade: 16S + KO\ncopy numbers", INK2, 8.5)
    rows = [
        ("s16fun\n(ours)", "#2a78d6", ["vsearch:\n25 most similar\nreference 16S", "evidence per KO:\nnearest 1 / 5 / 25 genomes,\nprevalence in genus..phylum,\nidentity, KO conservation",
                                       "boosted stacker\n-> probability\n+ expected copies", "expected-F1\nKO calls"]),
        ("PICRUSt2\n2.6.4", "#d55181", ["hmmalign + EPA-ng:\nplace in the\n16S tree", "castor hidden-state\nprediction: mp, emp_prob,\npic, subtree_average",
                                        "KO copy\nnumber", "calls: copies >= 1\n(or dev-tuned\nthreshold)"]),
        ("k-mer\nneural net", "#008300", ["1-5-mer\nfrequencies\n(1,364 features)", "2 x 512 ReLU network,\none sigmoid per KO\n(embeRNA design)", "probability", "calls >= 0.5"]),
        ("nearest-\nneighbour\nbaselines", "#eda100", ["vsearch:\nmost similar\nreference 16S", "copy nearest genome,\naverage 5 nearest,\nor genus core genome",
                                                      "Piphillin / Tax4Fun2 /\nMicFunPred-like:\nonly if >= 97 %", "calls"]),
    ]
    ys = [5.35, 3.85, 2.35, 0.85]
    widths = [1.75, 2.45, 1.85, 1.6]
    for (name, col, steps), y in zip(rows, ys):
        box(3.0, y, 1.25, 1.2, name, col, 9, True)
        x = 4.5
        for i, (st, w) in enumerate(zip(steps, widths)):
            box(x, y, w, 1.2, st, col, 8)
            arrow(x - 0.23, y + 0.6, x - 0.02, y + 0.6, col)
            x += w + 0.25
        arrow(2.5, 5.95, 2.98, y + 0.6)
        arrow(2.5, 2.1, 2.98, y + 0.6, INK2, ls="--")
        if name.startswith("s16fun") or name.startswith("PICRUSt2"):
            arrow(x - 0.23, y + 0.6, x + 0.02, y + 0.6, col)
            box(x + 0.05, y, 1.95, 1.2, ("community profile:\nsum of cells x\nexpected copies" if name.startswith("s16fun")
                                         else "community profile:\nsum of cells x\npredicted copies"), col, 8)
    ax.text(8.1, 0.2, "Held out per fold: whole species, genera, families, orders or classes (5-fold).   Dev families: method "
            "design and tuning.   Confirm families: scored once, with frozen settings.", ha="center", fontsize=9, color=INK2)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_primary(summary, prim, path):
    """Left: mean KO F1 per method and holdout level (full length). Middle and right: the co-primary endpoints,
    s16fun minus PICRUSt2 (best setting) with 95 % CIs (genera resampled), against their pre-registered bounds."""
    fig, axes = plt.subplots(1, 3, figsize=(17, 5), gridspec_kw={"width_ratios": [1.5, 1, 1]})
    a = axes[0]
    s = summary[summary.arm == "full"]
    x = np.arange(len(LEVELS))
    shown = [m for m in METHODS if m in set(s.method)]
    for i, m in enumerate(shown):
        lab, col, mk = METHODS[m]
        v = s[s.method == m].set_index("level").reindex(LEVELS).ko_f1
        a.plot(x + (i - len(shown) / 2) * 0.06, v, marker=mk, color=col, lw=1.5 if m in ("s16fun", "picrust2_best") else 0.9,
               ms=7 if m in ("s16fun", "picrust2_best") else 5, ls="--" if m == "null" else "-", label=lab)
    if "oracle5" in set(s.method):
        v = s[s.method == "oracle5"].set_index("level").reindex(LEVELS).ko_f1
        a.plot(x, v, color=INK, lw=1, ls=":", label="ceiling: 5 functionally nearest genomes (oracle, uses the truth)")
    a.set_xticks(x, [f"{lv}\nheld out" if lv != "random" else "random\n(species)" for lv in LEVELS])
    a.set_ylabel("KO F1 per genome (mean)", color=INK2)
    a.set_title("Accuracy as test genomes get more novel", color=INK, fontsize=11, loc="left")
    style(a)
    a.legend(fontsize=7.5, frameon=False, loc="lower left")
    if not prim.empty:
        y = np.arange(len(prim))
        for ax, col, bound, txt, title, ok in [
                (axes[1], "auroc", 0.05, " margin +0.05", "(a) per-KO AUROC: superiority", "auroc_superior"),
                (axes[2], "f1", -0.01, " bound -0.01", "(b) KO F1: non-inferiority", "f1_noninferior")]:
            d, lo, hi = prim[f"{col}_diff"], prim[f"{col}_ci_lo"], prim[f"{col}_ci_hi"]
            ax.axvline(0, color=INK2, lw=1)
            ax.axvline(bound, color="#d55181", lw=1, ls=":")
            ax.text(bound, len(prim) - 0.45, txt, color="#d55181", fontsize=8, va="bottom")
            ax.errorbar(d, y, xerr=[d - lo, hi - d], fmt="o", color="#2a78d6", ms=8, capsize=4, lw=1.6)
            for yi, dv, h, k in zip(y, d, hi, prim[ok]):
                ax.text(h + 0.003, yi, f"{dv:+.3f}  {'met' if k else 'not met'}", va="center", fontsize=9, color=INK)
            ax.set_yticks(y, [f"{lv} held out" for lv in prim.level])
            ax.invert_yaxis()
            ax.set_xlabel("s16fun minus PICRUSt2 (best setting)", color=INK2)
            ax.set_title(title, color=INK, fontsize=11, loc="left")
            style(ax)
            ax.grid(axis="x", color=GRID)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_novelty(novt, path):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
    for ax, by, xl in zip(axes, ["identity_bin", "nsti_bin"],
                          ["16S identity to the nearest training genome", "PICRUSt2 NSTI (tree distance to nearest reference)"]):
        d = novt[novt.bin_by == by]
        bins = list(dict.fromkeys(d.bin))
        for m in ["s16fun", "picrust2_best"]:
            lab, col, mk = METHODS[m]
            v = d[d.method == m].set_index("bin").reindex(bins)
            ax.plot(range(len(bins)), v.ko_f1, marker=mk, color=col, lw=1.8, ms=7, label=lab)
        n = d[d.method == "s16fun"].set_index("bin").reindex(bins).n
        ax.set_xticks(range(len(bins)), [f"{b}\nn={int(k) if k == k else 0}" for b, k in zip(bins, n)], fontsize=8)
        ax.set_xlabel(xl, color=INK2)
        style(ax)
    axes[0].set_ylabel("KO F1 per genome (mean)", color=INK2)
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle("Where each method wins: accuracy by novelty (full length, all holdout levels pooled)", color=INK, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_phyla(pg, genomes, path, min_n=30):
    d = pg[(pg.arm == "full") & (pg.level == "genus") & pg.method.isin(["s16fun", "picrust2_best"])]
    w = d.pivot_table(index="accession", columns="method", values="ko_f1").dropna()
    if w.empty:
        return
    w["phylum"] = genomes.loc[w.index, "phylum"]
    g = w.groupby("phylum").agg(n=("s16fun", "size"), s=("s16fun", "mean"), p=("picrust2_best", "mean"))
    g = g[g.n >= min_n].assign(diff=lambda x: x.s - x.p).sort_values("diff")
    fig, ax = plt.subplots(figsize=(9, max(3, 0.32 * len(g) + 1)))
    y = np.arange(len(g))
    ax.hlines(y, g.p, g.s, color=GRID, lw=3)
    ax.scatter(g.p, y, color=METHODS["picrust2_best"][1], marker="s", s=40, label=METHODS["picrust2_best"][0], zorder=3)
    ax.scatter(g.s, y, color=METHODS["s16fun"][1], marker="o", s=40, label=METHODS["s16fun"][0], zorder=3)
    ax.set_yticks(y, [f"{p} (n={n})" for p, n in zip(g.index, g.n)], fontsize=8)
    ax.set_xlabel("KO F1 per genome (mean), genus held out, full length", color=INK2)
    style(ax)
    ax.grid(axis="x", color=GRID)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_per_ko(A, path):
    fig, ax = plt.subplots(figsize=(6, 6))
    d = A[["picrust2_best", "s16fun"]].dropna()
    ax.scatter(d.picrust2_best, d.s16fun, s=3, alpha=0.3, color=METHODS["s16fun"][1], lw=0)
    ax.plot([0, 1], [0, 1], color=INK2, lw=1)
    ax.set_xlabel("per-KO AUROC, PICRUSt2 (best setting)", color=INK2)
    ax.set_ylabel("per-KO AUROC, s16fun", color=INK2)
    ax.text(0.02, 0.97, f"s16fun higher on {(d.s16fun > d.picrust2_best).mean():.0%} of {len(d):,} KOs", transform=ax.transAxes,
            va="top", fontsize=9, color=INK)
    style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_arms(summary, path):
    arms = [a for a in ["full", "V3V4", "V4"] if a in set(summary.arm)] + sorted(set(summary.arm) - {"full", "V3V4", "V4"})
    fig, axes = plt.subplots(1, len(arms), figsize=(4.5 * len(arms), 4.2), sharey=True, squeeze=False)
    for ax, arm in zip(axes[0], arms):
        s = summary[summary.arm == arm]
        for m in ["s16fun", "picrust2_best", "kmer_mlp", "knn5"]:
            if m not in set(s.method):
                continue
            lab, col, mk = METHODS[m]
            v = s[s.method == m].set_index("level").reindex(LEVELS).ko_f1
            ax.plot(range(len(LEVELS)), v, marker=mk, color=col, label=lab)
        ax.set_xticks(range(len(LEVELS)), LEVELS, fontsize=8)
        ax.set_title(f"{arm}", color=INK, fontsize=10, loc="left")
        style(ax)
    axes[0][0].set_ylabel("KO F1 per genome (mean)", color=INK2)
    axes[0][0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_ordination(ordi, path):
    spaces = ordi["spaces"]
    names = spaces[0]["names"][1:]
    fig, axes = plt.subplots(len(spaces), len(names), figsize=(3.9 * len(names), 3.9 * len(spaces)),
                             constrained_layout=True, squeeze=False)
    for r, sp in enumerate(spaces):
        n, C = sp["n"], sp["coords"]
        T = C[:n]
        for c, m in enumerate(names):
            ax = axes[r, c]
            P = C[(c + 1) * n:(c + 2) * n]
            ax.scatter(T[:, 0], T[:, 1], s=5, color=TRUTH, lw=0, label="truth")
            ax.scatter(P[:, 0], P[:, 1], s=5, color=METHODS[m][1], lw=0, alpha=0.7, label=METHODS[m][0])
            for i in sp["seg"]:
                ax.plot([T[i, 0], P[i, 0]], [T[i, 1], P[i, 1]], color=METHODS[m][1], lw=0.3, alpha=0.4)
            st = sp["stats"].set_index("method").loc[m]
            ax.text(0.02, 0.98, f"Mantel r {st.mantel_r:.2f}\nProcrustes m² {st.procrustes_m2:.2f}\n"
                    f"PERMANOVA R² {st.permanova_R2:.4f}\nANOSIM R {st.anosim_R:.4f}\n"
                    f"spread {st.spread_method / st.spread_truth:.2f}× truth", transform=ax.transAxes, va="top",
                    fontsize=7.5, color=INK, bbox=dict(fc="white", ec=GRID, alpha=0.85))
            ax.set_title(METHODS[m][0], fontsize=9, color=INK, loc="left")
            ax.set_xlabel(f"PCo1 ({sp['expl'][0]:.0%})", fontsize=8, color=INK2)
            ax.set_ylabel(f"PCo2 ({sp['expl'][1]:.0%})", fontsize=8, color=INK2)
            style(ax)
        axes[r, 0].annotate(sp["label"], xy=(0, 1.12), xycoords="axes fraction", fontsize=10, color=INK, fontweight="bold")
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_resources(res, path):
    keep = res[res.process.str.contains("BENCH_PREDICT|PICRUST2_FOLD|BENCH_MLP")].copy()
    if keep.empty:
        return
    lab = {"BENCH_PREDICT": "s16fun + nearest-neighbour baselines", "PICRUST2_FOLD": "PICRUSt2 (4 HSP settings)",
           "BENCH_MLP": "k-mer neural net"}
    keep["label"] = [next((v for k, v in lab.items() if k in p), p) for p in keep.process]
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.2))
    for ax, col, xl in zip(axes, ["median_task_min", "peak_gb"], ["median minutes per task", "peak memory per task (GB)"]):
        ax.barh(keep.label, keep[col], color=[METHODS["picrust2_best"][1] if "PICRUSt2" in x else METHODS["s16fun"][1]
                                              for x in keep.label])
        ax.set_xlabel(xl, color=INK2)
        style(ax)
        ax.grid(axis="x", color=GRID)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)

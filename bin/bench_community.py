#!/usr/bin/env python3
"""Simulated communities from held-out genomes: community KO profiles and differential abundance.

Usage: bench_community.py --data bench_data/ --tasks task_dir [...] [--level genus --methods mp,emp_prob,pic,subtree_average]
                          [--replicates 20 --samples 100 --pool 300 --responders 30 --fold-change 4 --seed 42]
                          --outdir community/

Communities are built only from genomes held out at --level (full-length 16S), so no method has them in its
reference. Each replicate draws a habitat pool of --pool genomes; each of --samples samples holds a random half
of the pool with log-normal cell abundances; in group B (half the samples) --responders pool genomes are
--fold-change times more abundant. Reads per genome = cells x true 16S copies. Each method turns reads back into
cells with its own predicted 16S copy number and into KO abundance with its own predicted KO copies:
  s16fun          expected copies, probability x expected copy number given presence (the method as shipped)
  s16fun_hard     its KO calls x expected copy number given presence (ablation: mean vs hard call)
  picrust2_<m>    PICRUSt2's predicted copy numbers, every HSP setting
Truth = true cells x true KO copies (PICRUSt2's table, capped at 10 as shipped).

Metrics per replicate (Sun et al. 2020 Microbiome 8:46): median per-sample Bray-Curtis and Spearman between
predicted and true relative KO abundance, plus the same with every method given the TRUE 16S copy numbers
(isolates the KO-content error from the 16S copy-number error); differential abundance A vs B: Spearman of
predicted and true per-KO log2 fold changes, and BH (q < 0.05) Wilcoxon calls scored against the truth's.

Outputs: community_reps.tsv (replicate x method), community_summary.tsv (mean and 95 % interval over replicates),
fig_community.png.
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import mannwhitneyu, spearmanr  # noqa: E402

from bench_figures import INK2, METHODS, style  # noqa: E402
from benchlib import load_data, load_pred  # noqa: E402


def collect(tasks, level, stems):
    """Predicted copies and 16S copies per held-out genome, for every method stem."""
    copies, ssu = {m: {} for m in stems}, {m: {} for m in stems}
    for t in tasks:
        info = pd.read_csv(t / "task.tsv", sep="\t")
        if info.empty or info.level.iloc[0] != level or info.arm.iloc[0] != "full":
            continue
        nov = pd.read_csv(t / "novelty.tsv", sep="\t", index_col=0).ssu_copies
        nsti = pd.read_csv(t / "nsti.tsv", sep="\t", index_col=0).ssu_copies
        for m, stem in stems.items():
            p = load_pred(t / f"pred_{stem}.npz")
            for i, q in enumerate(p["queries"]):
                copies[m][q] = p["copies"][i]
                ssu[m][q] = (nsti if stem.startswith("picrust2") else nov).get(q, 1.0)
    return copies, ssu


def bh(p):
    p = np.asarray(p)
    order = np.argsort(p)
    q = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1)
    return out


def da_calls(X, group):
    """KO x sample relative abundance -> BH-significant KOs (Wilcoxon rank-sum, two-sided)."""
    p = np.array([mannwhitneyu(x[group], x[~group]).pvalue if x.std() > 0 else 1.0 for x in X])
    return bh(p) < 0.05


def lfc(X, group):
    eps = 1e-9
    return np.log2(X[group].mean(0) + eps) - np.log2(X[~group].mean(0) + eps)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True)
    ap.add_argument("--tasks", nargs="+", required=True)
    ap.add_argument("--level", default="genus")
    ap.add_argument("--methods", default="mp,emp_prob,pic,subtree_average", help="PICRUSt2 HSP settings run")
    ap.add_argument("--replicates", type=int, default=20)
    ap.add_argument("--samples", type=int, default=100)
    ap.add_argument("--pool", type=int, default=300)
    ap.add_argument("--responders", type=int, default=30)
    ap.add_argument("--fold-change", type=float, default=4.0)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    genomes, kos, true_copies, _ = load_data(a.data)
    tasks = [Path(t) for t in a.tasks]
    stems = {"s16fun": "s16fun", **{f"picrust2_{m}": f"picrust2_{m}" for m in a.methods.split(",")}}
    copies, ssu = collect(tasks, a.level, stems)
    pool_all = sorted(set.intersection(*(set(c) for c in copies.values())))
    row = {g: i for i, g in enumerate(genomes.index)}
    T = np.minimum(true_copies[[row[g] for g in pool_all]].toarray(), 10).astype(np.float32)
    T16 = genomes.loc[pool_all, "ssu_copies"].to_numpy(float)
    pred = {m: np.vstack([copies[m][g] for g in pool_all]) for m in stems}
    calls = {}
    for t in tasks:   # s16fun hard calls, for the ablation
        info = pd.read_csv(t / "task.tsv", sep="\t")
        if not info.empty and info.level.iloc[0] == a.level and info.arm.iloc[0] == "full":
            p = load_pred(t / "pred_s16fun.npz")
            calls.update(zip(p["queries"], p["call"]))
    s16_hard = np.vstack([calls[g] for g in pool_all]) * expected_given_present(pred["s16fun"], tasks, a.level, pool_all)
    pred["s16fun_hard"] = s16_hard
    ssu["s16fun_hard"] = ssu["s16fun"]
    S16 = {m: np.maximum(np.array([ssu[m][g] for g in pool_all], float), 1) for m in pred}

    def profile(cells, C):
        X = cells @ C
        return X / np.maximum(X.sum(1, keepdims=True), 1e-12)

    rows = []
    for rep in range(a.replicates):
        pool = rng.choice(len(pool_all), min(a.pool, len(pool_all)), replace=False)
        resp = np.zeros(len(pool_all), bool)
        resp[rng.choice(pool, min(a.responders, len(pool)), replace=False)] = True
        group = np.arange(a.samples) >= a.samples // 2
        cells = np.zeros((a.samples, len(pool_all)))
        for s in range(a.samples):
            present = pool[rng.random(len(pool)) < 0.5]
            cells[s, present] = rng.lognormal(0, 1.5, len(present))
        cells[np.ix_(group, resp)] *= a.fold_change
        reads = cells * T16
        truth = profile(cells, T)
        keep = truth.std(0) > 0
        true_da = da_calls(truth[:, keep].T, group)
        true_lfc = lfc(truth[:, keep], group)
        for m, C in pred.items():
            for copies_from, est in (("own 16S copies", reads / S16[m]), ("true 16S copies", cells)):
                X = profile(est, C)
                bc = np.abs(X - truth).sum(1) / np.maximum((X + truth).sum(1), 1e-12)
                r = dict(replicate=rep, method=m, copies_from=copies_from, bray_curtis=np.median(bc),
                         spearman=np.median([spearmanr(X[i], truth[i]).statistic for i in range(a.samples)]))
                if copies_from == "own 16S copies":
                    pd_ = da_calls(X[:, keep].T, group)
                    tp, fp, fn = (pd_ & true_da).sum(), (pd_ & ~true_da).sum(), (~pd_ & true_da).sum()
                    r.update(lfc_spearman=spearmanr(lfc(X[:, keep], group), true_lfc).statistic,
                             da_true=int(true_da.sum()), da_predicted=int(pd_.sum()), da_precision=tp / max(tp + fp, 1),
                             da_recall=tp / max(tp + fn, 1), da_f1=2 * tp / max(2 * tp + fp + fn, 1))
                rows.append(r)
    reps = pd.DataFrame(rows)
    reps.to_csv(out / "community_reps.tsv", sep="\t", index=False)
    cols = ["bray_curtis", "spearman", "lfc_spearman", "da_precision", "da_recall", "da_f1"]
    summ = reps.groupby(["copies_from", "method"])[cols].agg(["mean", lambda x: np.percentile(x, 2.5),
                                                               lambda x: np.percentile(x, 97.5)])
    summ.columns = [f"{c}_{s}" for c, s in zip(np.repeat(cols, 3), ["mean", "lo", "hi"] * len(cols))]
    summ.reset_index().to_csv(out / "community_summary.tsv", sep="\t", index=False)

    own = reps[reps.copies_from == "own 16S copies"]
    order = own.groupby("method").bray_curtis.mean().sort_values().index.tolist()
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    for ax, col, lab in zip(axes, ["bray_curtis", "spearman", "lfc_spearman"],
                            ["Bray-Curtis to the true profile (lower is better)", "Spearman with the true profile",
                             "Spearman of log2 fold changes, B vs A"]):
        for i, m in enumerate(order):
            v = own[own.method == m][col]
            colr = METHODS["s16fun"][1] if m.startswith("s16fun") else METHODS["picrust2_best"][1]
            ax.boxplot(v, positions=[i], widths=0.6, patch_artist=True, showfliers=False,
                       boxprops=dict(facecolor=colr + ("" if m == "s16fun" else "66"), color=INK2), medianprops=dict(color=INK2))
        ax.set_xticks(range(len(order)), [m.replace("picrust2_", "PICRUSt2 ").replace("_", " ") for m in order],
                      rotation=35, ha="right", fontsize=8)
        ax.set_ylabel(lab, color=INK2, fontsize=9)
        style(ax)
    fig.tight_layout()
    fig.savefig(out / "fig_community.png", dpi=150)
    print(f"{len(pool_all)} held-out genomes; {a.replicates} replicates")


def expected_given_present(expected, tasks, level, pool_all):
    """Expected copies given presence = expected copies / probability (where the probability is not ~0)."""
    prob = {}
    for t in tasks:
        info = pd.read_csv(t / "task.tsv", sep="\t")
        if not info.empty and info.level.iloc[0] == level and info.arm.iloc[0] == "full":
            p = load_pred(t / "pred_s16fun.npz")
            prob.update(zip(p["queries"], p["score"]))
    Pm = np.vstack([prob[g] for g in pool_all])
    return np.where(Pm > 0.01, np.minimum(expected / np.maximum(Pm, 1e-6), 10), 1.0)


if __name__ == "__main__":
    main()

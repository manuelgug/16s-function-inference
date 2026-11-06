#!/usr/bin/env python3
"""Simulated communities from held-out genomes: community KO profiles and differential abundance.

Usage: bench_community.py --data bench_data/ --tasks task_dir [...] [--level genus --arms V4,full
                          --methods mp,emp_prob,pic,subtree_average --replicates 10 --seed 42] --outdir community/

Communities (bin/community.py) are built only from genomes held out at --level, so no method has them in its
reference; one ASV per genome = its sequence in the task's arm (V4 amplicon: what users have). Habitats: gut,
marine, soil designs. Each method turns reads into cells with its own predicted 16S copy number and into KO
abundance with its own predicted KO copies:
  s16fun          expected copies (probability x expected copy number when present), learned 16S copies
  s16fun_hard     its KO calls x expected copy number when present (ablation: expectation vs hard call)
  picrust2_<m>    PICRUSt2's predicted copy numbers and 16S copies, every HSP setting
Each method is also scored with the TRUE 16S copy numbers (isolates the KO-content error).
Truth = true cells x true KO copies (PICRUSt2's table, capped at 10 as shipped).

Outputs: community_reps.tsv (arm x habitat x replicate x method x 16S source), community_summary.tsv,
fig_community.png.
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from bench_figures import INK2, METHODS, style  # noqa: E402
from benchlib import load_data, load_pred  # noqa: E402
from community import HABITATS, profile, score, simulate  # noqa: E402


def collect(tasks, level, arm, stems):
    """Per held-out genome of one level x arm: predicted copies and 16S copies per method, s16fun probabilities/calls."""
    out = {m: {} for m in stems}
    ssu = {m: {} for m in stems}
    prob, call = {}, {}
    for t in tasks:
        info = pd.read_csv(t / "task.tsv", sep="\t")
        if info.empty or info.level.iloc[0] != level or info.arm.iloc[0] != arm:
            continue
        nov = pd.read_csv(t / "novelty.tsv", sep="\t", index_col=0).ssu_copies
        nsti = pd.read_csv(t / "nsti.tsv", sep="\t", index_col=0).ssu_copies
        for m, stem in stems.items():
            p = load_pred(t / f"pred_{stem}.npz")
            for i, q in enumerate(p["queries"]):
                out[m][q] = p["copies"][i]
                ssu[m][q] = (nsti if stem.startswith("picrust2") else nov).get(q, 1.0)
            if stem == "s16fun":
                prob.update(zip(p["queries"], p["score"]))
                call.update(zip(p["queries"], p["call"]))
    return out, ssu, prob, call


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True)
    ap.add_argument("--tasks", nargs="+", required=True)
    ap.add_argument("--level", default="genus")
    ap.add_argument("--arms", default="V4,full")
    ap.add_argument("--methods", default="mp,emp_prob,pic,subtree_average", help="PICRUSt2 HSP settings run")
    ap.add_argument("--replicates", type=int, default=10)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    genomes, kos, true_copies, _ = load_data(a.data)
    row = {g: i for i, g in enumerate(genomes.index)}
    tasks = [Path(t) for t in a.tasks]
    stems = {"s16fun": "s16fun", **{f"picrust2_{m}": f"picrust2_{m}" for m in a.methods.split(",")}}

    rows = []
    for arm in a.arms.split(","):
        copies, ssu, prob, call = collect(tasks, a.level, arm, stems)
        if not prob:
            continue
        pool = sorted(set.intersection(*(set(c) for c in copies.values())))
        T = np.minimum(true_copies[[row[g] for g in pool]].toarray(), 10).astype(np.float32)
        T16 = genomes.loc[pool, "ssu_copies"].to_numpy(float)
        pred = {m: np.vstack([copies[m][g] for g in pool]) for m in stems}
        P = np.vstack([prob[g] for g in pool])
        cond = np.where(P > 0.01, np.minimum(pred["s16fun"] / np.maximum(P, 1e-6), 10), 1.0)
        pred["s16fun_hard"] = np.vstack([call[g] for g in pool]) * cond
        ssu["s16fun_hard"] = ssu["s16fun"]
        S16 = {m: np.maximum(np.array([ssu[m][g] for g in pool], float), 1) for m in pred}
        for habitat in HABITATS:
            for rep in range(a.replicates):
                cells, group, _ = simulate(len(pool), habitat, rng)
                reads = cells * T16
                truth = profile(cells, T)
                for m, C in pred.items():
                    for src, est in (("own 16S copies", reads / S16[m]), ("true 16S copies", cells)):
                        r = score(profile(est, C), truth, group)
                        rows.append(dict(arm=arm, habitat=habitat, replicate=rep, method=m, copies_from=src, **r))
    reps = pd.DataFrame(rows)
    reps.to_csv(out / "community_reps.tsv", sep="\t", index=False)
    cols = ["bray_curtis", "spearman", "lfc_spearman", "da_precision", "da_recall", "da_f1"]
    g = reps.groupby(["arm", "habitat", "copies_from", "method"])[cols]
    summ = g.mean().add_suffix("_mean").join(g.quantile(0.025).add_suffix("_lo")).join(g.quantile(0.975).add_suffix("_hi"))
    summ.reset_index().to_csv(out / "community_summary.tsv", sep="\t", index=False)

    own = reps[reps.copies_from == "own 16S copies"]
    arm0 = a.arms.split(",")[0] if a.arms.split(",")[0] in set(own.arm) else own.arm.iloc[0]
    own = own[own.arm == arm0]
    order = own.groupby("method").bray_curtis.mean().sort_values().index.tolist()
    fig, axes = plt.subplots(1, len(HABITATS), figsize=(5 * len(HABITATS), 4.4), sharey=True)
    for ax, habitat in zip(axes, HABITATS):
        for i, m in enumerate(order):
            v = own[(own.method == m) & (own.habitat == habitat)].bray_curtis
            colr = METHODS["s16fun"][1] if m.startswith("s16fun") else METHODS["picrust2_best"][1]
            ax.boxplot(v, positions=[i], widths=0.6, patch_artist=True, showfliers=False,
                       boxprops=dict(facecolor=colr if m == "s16fun" else colr + "55", color=INK2),
                       medianprops=dict(color=INK2))
        ax.set_xticks(range(len(order)), [m.replace("picrust2_", "PICRUSt2 ").replace("_", " ") for m in order],
                      rotation=35, ha="right", fontsize=8)
        ax.set_title(f"{habitat} ({HABITATS[habitat]['richness']} ASVs per sample)", fontsize=10, loc="left")
        style(ax)
    axes[0].set_ylabel(f"Bray-Curtis to the true KO profile ({arm0} ASVs; lower is better)", color=INK2, fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "fig_community.png", dpi=150)
    print(f"{len(reps)} rows")


if __name__ == "__main__":
    main()

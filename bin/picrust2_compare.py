#!/usr/bin/env python3
"""Head-to-head: our 16S models against PICRUSt2 (M2), per genome, on the same genomes and targets.

Usage: picrust2_compare.py --pred combined_KO_predicted.tsv.gz --nsti combined_marker_predicted_and_nsti.tsv.gz
           --refs reference_genomes.txt --ko-names ko_info.tsv.gz --targets targets/ --aim1 full_random/
           --genomes genomes.tsv [--cpus 16] [--seed 42] [--boot 1000] [--max-ordination 1500]
           --outdir picrust2_comparison

Our predictions are the out-of-fold predictions of the full-length, random-split Aim 1 run: every
genome is predicted by models that never saw it. PICRUSt2 2.6.4 runs with its default reference.
The fair comparison set is the genomes NOT in PICRUSt2's reference (otherwise PICRUSt2 looks up
the genome's own annotation). KOs: the Aim 1 KO set (prevalence filter); PICRUSt2 calls a KO
present at predicted copy number >= 1, our models at probability >= 0.5. Modules: completeness from
predicted KOs (PICRUSt2: kegg-pathways-completeness on its full predicted KO set, as for the truth).
Outputs: per_genome.tsv, paired.tsv, by_novelty.tsv, by_taxon.tsv, per_ko.tsv, ordination_stats.tsv,
ko_calls.tsv.gz, module_calls.tsv.gz, figures and picrust2_comparison.md (+ .html via render_html.py).
"""
import argparse
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.linalg import eigh  # noqa: E402
from scipy.spatial import procrustes  # noqa: E402
from scipy.spatial.distance import cdist  # noqa: E402
from scipy.stats import rankdata, wilcoxon  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

from aim1_report import GRID, ID_BINS, ID_LABELS, INK2, boot_ci, md_table, quote, style  # noqa: E402
from build_targets import module_matrix  # noqa: E402
from picrust2_eval import in_reference  # noqa: E402

# PICRUSt2 takes dataviz palette slot 5 (magenta); our models keep their Aim 1 colours.
METHODS = {"M2_picrust2": ("PICRUSt2", "#d55181"), "M1_knn1": ("1-NN", "#2a78d6"),
           "M1_knn5": ("5-NN", "#eb6834"), "M3_rf": ("random forest", "#1baf7a"), "M4_lgbm": ("LightGBM", "#eda100")}
OURS = ["M1_knn1", "M1_knn5", "M3_rf", "M4_lgbm"]
TRUTH_COLOR = "#b9b8b3"
LOWER_BETTER = {"module_mae"}


def ko_scores(Y, P):
    tp, fp, fn = (P & Y).sum(1), (P & ~Y).sum(1), (~P & Y).sum(1)
    return {"ko_f1": 2 * tp / np.maximum(2 * tp + fp + fn, 1), "ko_jaccard": tp / np.maximum(tp + fp + fn, 1),
            "ko_precision": tp / np.maximum(tp + fp, 1), "ko_recall": tp / np.maximum(tp + fn, 1),
            "ko_n_predicted": P.sum(1), "ko_tp": tp, "ko_fp": fp, "ko_fn": fn}


def jaccard_dist(A):
    A = A.astype("float32")
    inter = A @ A.T
    s = A.sum(1)
    union = s[:, None] + s[None, :] - inter
    return np.where(union > 0, 1 - inter / np.maximum(union, 1), 0.0).astype("float64")


def bray_curtis_dist(M):
    s = M.sum(1)
    return cdist(M, M, "cityblock") / np.maximum(s[:, None] + s[None, :], 1e-12)


def pcoa(D, k=2):
    """Classical scaling: top-k coordinates and their share of total inertia (trace of the centred matrix)."""
    D2 = D ** 2
    B = -0.5 * (D2 - D2.mean(0)[None, :] - D2.mean(1)[:, None] + D2.mean())
    n = len(B)
    w, v = eigh(B, subset_by_index=[n - k, n - 1])
    w, v = w[::-1], v[:, ::-1]
    return v * np.sqrt(np.maximum(w, 0)), w / np.trace(B)


def group_tests(D, n, rng, perms=999):
    """PERMANOVA, ANOSIM and PERMDISP for two groups of n paired profiles (rows 0..n-1 vs n..2n-1, row i
    paired with n+i: same genome). Permutations swap the two labels within each genome (strata = genome).
    Returns effect sizes, permutation p-values and each group's mean distance to its centroid."""
    N = 2 * n
    Q = D ** 2
    r = rankdata(D[np.triu_indices(N, 1)])
    R = np.zeros_like(D)
    R[np.triu_indices(N, 1)] = r
    R = R + R.T
    swap = rng.random((perms, n)) < 0.5
    A = np.zeros((N, perms + 1))
    A[n:, 0] = 1  # observed: second block is group 1
    A[:n, 1:] = swap.T
    A[n:, 1:] = ~swap.T
    n1 = n0 = n
    ones_q, ones_r = Q.sum(1), R.sum(1)

    def within(M, rowsum):
        MA = M @ A
        aMa = (A * MA).sum(0)
        bMb = rowsum.sum() - 2 * (rowsum @ A) + aMa
        return MA, aMa, bMb
    QA, aQa, bQb = within(Q, ones_q)
    ss_t = Q.sum() / (2 * N)
    ss_w = aQa / (2 * n1) + bQb / (2 * n0)
    F = (ss_t - ss_w) / (ss_w / (N - 2))
    _, aRa, bRb = within(R, ones_r)
    w_pairs = n1 * (n1 - 1) / 2 + n0 * (n0 - 1) / 2
    rw = (aRa + bRb) / 2 / w_pairs
    rb = (r.sum() - (aRa + bRb) / 2) / (n1 * n0)
    Ran = (rb - rw) / (len(r) / 2)
    QB = ones_q[:, None] - QA  # Q @ (1 - a)
    z2 = A * (QA / n1 - aQa / (2 * n1 ** 2)) + (1 - A) * (QB / n0 - bQb / (2 * n0 ** 2))
    z = np.sqrt(np.maximum(z2, 0))
    m1, m0 = (A * z).sum(0) / n1, ((1 - A) * z).sum(0) / n0
    zbar = z.mean(0)
    ssb = n1 * (m1 - zbar) ** 2 + n0 * (m0 - zbar) ** 2
    ssw = (A * (z - m1) ** 2).sum(0) + ((1 - A) * (z - m0) ** 2).sum(0)
    Fd = ssb / (ssw / (N - 2))
    # ties count as >= (tolerance: identical groups give equal statistics up to rounding)
    pval = lambda x: (1 + (x[1:] >= x[0] - 1e-9 * max(1.0, abs(x[0]))).sum()) / (perms + 1)  # noqa: E731
    return {"permanova_R2": (ss_t - ss_w[0]) / ss_t, "permanova_p": pval(F), "anosim_R": Ran[0],
            "anosim_p": pval(Ran), "dispersion_group0": m0[0], "dispersion_group1": m1[0],
            "permdisp_p": pval(Fd)}


def paired_rows(pg, subset, metrics, rng, n_boot):
    rows = []
    for m in metrics:
        sign = -1 if m in LOWER_BETTER else 1
        for model in OURS:
            if f"{m}__{model}" not in pg:  # e.g. LightGBM is fitted for modules and COGs only
                continue
            a, b = pg.loc[subset, f"{m}__{model}"], pg.loc[subset, f"{m}__M2_picrust2"]
            ok = a.notna() & b.notna()
            if ok.sum() < 10:
                continue
            d = sign * (a[ok] - b[ok]).to_numpy()
            lo, hi = boot_ci(d, rng, n_boot)
            rows.append({"metric": m, "ours": METHODS[model][0], "n": int(ok.sum()), "ours_mean": a[ok].mean(),
                         "picrust2_mean": b[ok].mean(), "diff (ours better > 0)": d.mean(), "ci_lo": lo, "ci_hi": hi,
                         "wilcoxon_p": wilcoxon(d).pvalue if np.any(d) else 1.0,
                         "genomes_ours_better": (d > 0).mean(), "genomes_picrust2_better": (d < 0).mean()})
    return pd.DataFrame(rows)


def per_ko(Y, scores, calls, kos):
    """Per-KO AUROC (ranking) and F1 (calls) for every method, over the given genomes."""
    rows = []
    for j, ko in enumerate(kos):
        y = Y[:, j]
        if not 0 < y.sum() < len(y):
            continue
        r = {"ko": ko, "prevalence": y.mean()}
        for m in scores:
            s, c = scores[m][:, j], calls[m][:, j]
            r[f"auroc__{m}"] = roc_auc_score(y, s) if np.ptp(s) else 0.5
            tp, fp, fn = (c & y).sum(), (c & ~y).sum(), (~c & y).sum()
            r[f"f1__{m}"] = 2 * tp / max(2 * tp + fp + fn, 1)
        rows.append(r)
    return pd.DataFrame(rows)


def fig_novelty(pg, sub, out):
    """Accuracy by PICRUSt2 NSTI quintile and by our nearest-training 16S identity, unseen genomes."""
    d = pg.loc[sub].copy()
    d["nsti_bin"] = pd.qcut(d["nsti"], 5, duplicates="drop")
    d["id_bin"] = pd.cut(d["nn_identity"], ID_BINS, labels=ID_LABELS, right=False)
    panels = [("nsti_bin", "ko_f1", "KO set F1 (higher is better)", "PICRUSt2 NSTI (quintiles)"),
              ("id_bin", "ko_f1", "KO set F1 (higher is better)", "16S identity to nearest training genome (%)"),
              ("nsti_bin", "module_mae", "Module completeness MAE (lower is better)", "PICRUSt2 NSTI (quintiles)")]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True)
    rows = []
    for ax, (b, m, title, xl) in zip(axes, panels):
        for model, (label, color) in METHODS.items():
            col = f"{m}__{model}"
            if col not in d or d[col].isna().all():
                continue
            s = d.groupby(b, observed=False)[col].agg(["mean", "size"])
            x = np.arange(len(s))
            ax.plot(x, s["mean"], marker="o", color=color, label=label)
            rows += [{"bin_by": b, "bin": str(i), "metric": m, "method": label, "mean": v, "n": n}
                     for i, v, n in zip(s.index, s["mean"], s["size"])]
        cats = [f"{i.left:.2f}-{i.right:.2f}" if hasattr(i, "left") else str(i) for i in s.index]
        ax.set_xticks(x, [f"{c}\nn={n}" for c, n in zip(cats, s["size"])], fontsize=8)
        ax.set_title(title)
        ax.set_xlabel(xl)
    axes[0].set_ylabel("mean per-genome value")
    h, lab = axes[2].get_legend_handles_labels()
    fig.legend(h, lab, loc="outside upper center", ncol=len(lab))
    fig.savefig(out / "fig1_accuracy_by_novelty.png", dpi=150)
    plt.close(fig)
    return pd.DataFrame(rows)


def fig_per_ko(pk, out):
    """Left: ranking (AUROC) of the best-ranking model; right: calls (F1) of the best-calling model."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), constrained_layout=True)
    for ax, (score, model, lo) in zip(axes, [("auroc", "M3_rf", 0.2), ("f1", "M1_knn5", 0.0)]):
        x, y = pk[f"{score}__M2_picrust2"], pk[f"{score}__{model}"]
        ax.scatter(x, y, s=6, alpha=0.4, color=METHODS[model][1], lw=0)
        ax.plot([0, 1], [0, 1], color=INK2, lw=1)
        ax.set_xlim(lo, 1.01)
        ax.set_ylim(lo, 1.01)
        ax.set_aspect("equal")
        ax.set_xlabel(f"PICRUSt2 per-KO {score.upper()}")
        ax.set_ylabel(f"{METHODS[model][0]} per-KO {score.upper()}")
        ax.set_title(f"{score.upper()}: {METHODS[model][0]} better on {(y > x).mean():.0%}, PICRUSt2 on "
                     f"{(y < x).mean():.0%} of {len(pk)} KOs", fontsize=10)
    fig.savefig(out / "fig2_per_ko_auroc.png", dpi=150)
    plt.close(fig)


def ordination(truth, preds, dist, label, rng, max_n, seg_n=120):
    """Joint PCoA of truth + every method's profiles for the same genomes; per-method agreement with truth."""
    idx = np.arange(len(truth))
    if len(idx) > max_n:  # shortcut: dense n^2 PCoA; subsample when many genomes (full scale)
        idx = np.sort(rng.choice(idx, max_n, replace=False))
    names = ["truth"] + list(preds)
    stacked = np.vstack([truth[idx]] + [preds[m][idx] for m in preds])
    D = dist(stacked)
    coords, expl = pcoa(D, 2)
    n = len(idx)
    Dt = D[:n, :n]
    iu = np.triu_indices(n, 1)
    t10, _ = pcoa(Dt, 10)
    stats, coord_rows = [], []
    for k, m in enumerate(names):
        blk = slice(k * n, (k + 1) * n)
        coord_rows.append(pd.DataFrame({"profile": m, "genome_index": idx, "PCo1": coords[blk, 0],
                                        "PCo2": coords[blk, 1], "space": label}))
        if m == "truth":
            continue
        Dm = D[blk, blk]
        m10, _ = pcoa(Dm, 10)
        both = np.r_[np.arange(n), np.arange(k * n, (k + 1) * n)]
        gt = group_tests(D[np.ix_(both, both)], n, rng)
        stats.append({"space": label, "method": METHODS[m][0], "genomes": n,
                      "mantel_r (between-genome distances vs truth)": np.corrcoef(Dt[iu], Dm[iu])[0, 1],
                      "procrustes_m2 (0 = same configuration)": procrustes(t10, m10)[2],
                      "mean distance to own truth": np.diag(D[blk, :n]).mean(),
                      "PERMANOVA R2": gt["permanova_R2"], "PERMANOVA p": gt["permanova_p"],
                      "ANOSIM R": gt["anosim_R"], "ANOSIM p": gt["anosim_p"],
                      "spread truth": gt["dispersion_group0"], "spread method": gt["dispersion_group1"],
                      "PERMDISP p": gt["permdisp_p"]})
    # our models against PICRUSt2 directly (same genomes, paired the same way)
    pk_i = names.index("M2_picrust2") if "M2_picrust2" in names else None
    vs = []
    if pk_i is not None:
        for k, m in enumerate(names):
            if m in ("truth", "M2_picrust2"):
                continue
            both = np.r_[np.arange(pk_i * n, (pk_i + 1) * n), np.arange(k * n, (k + 1) * n)]
            gt = group_tests(D[np.ix_(both, both)], n, rng)
            vs.append({"space": label, "comparison": f"{METHODS[m][0]} vs PICRUSt2", "genomes": n,
                       "PERMANOVA R2": gt["permanova_R2"], "PERMANOVA p": gt["permanova_p"],
                       "ANOSIM R": gt["anosim_R"], "ANOSIM p": gt["anosim_p"],
                       "spread PICRUSt2": gt["dispersion_group0"], "spread ours": gt["dispersion_group1"],
                       "PERMDISP p": gt["permdisp_p"]})
    seg = rng.choice(n, min(seg_n, n), replace=False)
    return coords, expl, n, pd.DataFrame(stats), pd.concat(coord_rows), seg, pd.DataFrame(vs)


def fig_ordination(rows, out):
    fig, axes = plt.subplots(len(rows), 4, figsize=(15, 3.9 * len(rows)), constrained_layout=True,
                             squeeze=False)
    for r, (label, coords, expl, n, st, methods, seg) in enumerate(rows):
        T = coords[:n]
        for c, m in enumerate(methods):
            ax = axes[r, c]
            P = coords[(c + 1) * n:(c + 2) * n]
            ax.scatter(T[:, 0], T[:, 1], s=5, color=TRUTH_COLOR, lw=0, label="truth (eggNOG annotation)")
            ax.scatter(P[:, 0], P[:, 1], s=5, color=METHODS[m][1], lw=0, alpha=0.7, label=METHODS[m][0])
            for i in seg:
                ax.plot([T[i, 0], P[i, 0]], [T[i, 1], P[i, 1]], color=METHODS[m][1], lw=0.5, alpha=0.5)
            s = st[st.method == METHODS[m][0]].iloc[0]
            ax.set_title(f"{label}: {METHODS[m][0]}", fontsize=10)
            ax.text(0.02, 0.98, f"Mantel r {s['mantel_r (between-genome distances vs truth)']:.2f} · Procrustes m² "
                    f"{s['procrustes_m2 (0 = same configuration)']:.2f}\n"
                    f"PERMANOVA R² {s['PERMANOVA R2']:.4f} (p {s['PERMANOVA p']:.3f})\n"
                    f"ANOSIM R {s['ANOSIM R']:.4f} (p {s['ANOSIM p']:.3f})\n"
                    f"spread {s['spread method'] / s['spread truth']:.2f}× truth (PERMDISP p {s['PERMDISP p']:.3f})",
                    transform=ax.transAxes, va="top", ha="left", fontsize=7.5, color=INK2,
                    bbox=dict(boxstyle="round,pad=0.3", fc="#fcfcfb", ec=GRID, alpha=0.9))
            ax.set_xlabel(f"PCo1 ({expl[0]:.0%})")
            if c == 0:
                ax.set_ylabel(f"PCo2 ({expl[1]:.0%})")
            if r > 0 or c > 0:
                ax.sharex(axes[r, 0])
                ax.sharey(axes[r, 0])
    handles = [plt.Line2D([], [], marker="o", ls="", color=TRUTH_COLOR, label="truth (eggNOG annotation)")]
    handles += [plt.Line2D([], [], marker="o", ls="", color=METHODS[m][1], label=METHODS[m][0]) for m in rows[0][5]]
    fig.legend(handles=handles, loc="outside upper center", ncol=len(handles))
    fig.savefig(out / "fig3_ordination.png", dpi=150)
    plt.close(fig)


DIFFERENCES = [
    ("Question answered", "Function of one genome (or one ASV) from its 16S", "Same per sequence; built mainly to "
     "predict community (metagenome) profiles from ASV tables, with 16S copy-number normalisation"),
    ("Reference genomes", "GTDB R232 representatives passing CheckM2 >= 90 % / <= 5 % and a full-length 16S "
     "(this run: {n_ref} genomes; each genome predicted from the other ~80 %)",
     "PICRUSt2-SC: GTDB r214 representatives with a high-quality 16S: 26,868 bacteria + 1,002 archaea"),
    ("Function annotation of the reference", "eggNOG-mapper 2.1.15, eggNOG DB 5.0.2, pyrodigal genes",
     "eggNOG-mapper 2.1.12, eggNOG DB 5.0.2, Prodigal genes (the same pipeline as our ground truth)"),
    ("16S processing", "SINA alignment to SILVA 144; aligned positions become features",
     "HMMER alignment to the reference MSA; EPA-ng placement into the GTDB-based 16S tree"),
    ("Prediction", "k-NN on 16S identity (copy the closest genomes) or learned models (random forest, LightGBM) "
     "on alleles", "Hidden-state prediction along the tree (castor, maximum parsimony by default)"),
    ("Output per genome", "KO presence probability; KEGG module completeness; COG-category composition",
     "KO (and EC, ...) copy numbers; pathways via MinPath for communities"),
    ("Amplicons", "One model per primer pair (V4, V3-V4, full length), trained on the in-silico amplicon",
     "Any fragment is placed on the same full-length reference tree"),
    ("Novelty measure", "16S identity to the nearest training genome",
     "NSTI: tree distance to the nearest reference genome; sequences with NSTI > 2 dropped"),
    ("Validation built in", "Clade-blocked cross-validation (genus to class held out)", "None per run; NSTI as a warning"),
    ("Tunable to a new reference", "Yes: rerun on any genome set (e.g. all 55k GTDB R232 genomes)",
     "Yes, with a custom reference package (tree, HMM, trait tables)"),
]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for x in ["--pred", "--nsti", "--refs", "--ko-names", "--targets", "--aim1", "--genomes", "--outdir"]:
        ap.add_argument(x, required=True)
    ap.add_argument("--cpus", type=int, default=1)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--max-ordination", type=int, default=1500)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    style()
    pdir = Path(a.aim1) / "predictions"

    # ---- truth and predictions on one genome x target grid ----
    ours_ko = {m: pd.read_parquet(pdir / f"{m}_ko.parquet") for m in OURS if (pdir / f"{m}_ko.parquet").exists()}
    ours_mod = {m: pd.read_parquet(pdir / f"{m}_module.parquet") for m in OURS if (pdir / f"{m}_module.parquet").exists()}
    genomes = ours_ko["M1_knn5"].index
    kos = ours_ko["M1_knn5"].columns
    mods = ours_mod["M1_knn5"].columns
    t = Path(a.targets)
    Y = pd.read_parquet(t / "ko_presence.parquet").reindex(index=genomes, columns=kos).fillna(0).to_numpy() >= 0.5
    Ym = pd.read_parquet(t / "module_completeness.parquet").reindex(index=genomes, columns=mods).fillna(0).to_numpy()

    pc = pd.read_csv(a.pred, sep="\t", index_col=0)
    pc.columns = pc.columns.str.removeprefix("ko:")
    dropped = genomes.difference(pc.index)  # PICRUSt2 drops sequences with NSTI > 2 (no prediction = all absent)
    pc = pc.reindex(genomes).fillna(0)
    pc_kos = {g: set(pc.columns[row >= 1]) for g, row in zip(pc.index, pc.to_numpy())}
    with tempfile.TemporaryDirectory() as tmp:
        pc_mod = module_matrix(pc_kos, list(genomes), tmp, a.cpus).reindex(columns=mods).fillna(0.0)
    scores = {"M2_picrust2": pc.reindex(columns=kos).fillna(0).to_numpy()}
    scores |= {m: p.to_numpy() for m, p in ours_ko.items()}
    calls = {"M2_picrust2": scores["M2_picrust2"] >= 1} | {m: scores[m] >= 0.5 for m in ours_ko}
    mod_pred = {"M2_picrust2": pc_mod.to_numpy()} | {m: p.to_numpy() for m, p in ours_mod.items()}

    # ---- per genome ----
    meta = pd.read_csv(a.genomes, sep="\t").set_index("accession").reindex(genomes)
    nsti = pd.read_csv(a.nsti, sep="\t", index_col=0).reindex(genomes)
    near = pd.read_csv(Path(a.aim1) / "nearest.tsv", sep="\t").set_index("accession").reindex(genomes)
    pg = pd.DataFrame(index=genomes)
    pg["domain"], pg["phylum"], pg["species"] = meta["domain"], meta["phylum"], meta["species"]
    pg["in_picrust2_reference"] = in_reference(genomes, a.refs)
    pg["picrust2_dropped"] = genomes.isin(dropped)
    pg["nsti"] = nsti["metadata_NSTI"]
    pg["picrust2_closest_reference"] = nsti["closest_reference_genome"]
    pg["nn_identity"] = near["nn_identity"]
    pg["ko_n_true"] = Y.sum(1)
    for m in calls:
        for k, v in ko_scores(Y, calls[m]).items():
            pg[f"{k}__{m}"] = v
    for m in mod_pred:
        pg[f"module_mae__{m}"] = np.abs(mod_pred[m] - Ym).mean(1)
    pg.index.name = "accession"
    pg.to_csv(out / "per_genome.tsv", sep="\t")
    unseen = ~pg["in_picrust2_reference"]

    # ---- 1. head to head ----
    metrics = ["ko_f1", "ko_jaccard", "ko_precision", "ko_recall", "module_mae"]
    paired = pd.concat([paired_rows(pg, s, metrics, rng, a.boot).assign(genomes=name)
                        for name, s in [("not in PICRUSt2 reference", unseen),
                                        ("in PICRUSt2 reference (leak)", ~unseen)]])
    paired.to_csv(out / "paired.tsv", sep="\t", index=False)
    means = pd.DataFrame({METHODS[m][0]: [pg.loc[unseen, f"{k}__{m}"].mean() if f"{k}__{m}" in pg else np.nan
                                          for k in metrics + ["ko_n_predicted"]]
                          for m in METHODS}, index=metrics + ["ko_n_predicted"])
    means.loc["ko_n_true"] = pg.loc[unseen, "ko_n_true"].mean()
    means = means.reset_index().rename(columns={"index": "metric (mean over unseen genomes)"})

    f1 = paired[(paired.metric == "ko_f1") & (paired.genomes == "not in PICRUSt2 reference")].set_index("ours")
    best = f1["ours_mean"].idxmax()
    b = f1.loc[best]
    verdict = ("beats" if b.ci_lo > 0 else "is beaten by" if b.ci_hi < 0 else "is not distinguishable from")
    mm = paired[(paired.metric == "module_mae") & (paired.genomes == "not in PICRUSt2 reference")].set_index("ours")
    mbest = mm["ours_mean"].idxmin()
    mb = mm.loc[mbest]
    mverdict = ("beats" if mb.ci_lo > 0 else "is beaten by" if mb.ci_hi < 0 else "is not distinguishable from")

    # ---- 2. where ----
    nov = fig_novelty(pg, unseen, out)
    nov.to_csv(out / "by_novelty.tsv", sep="\t", index=False)
    tax = []
    for col in ["domain", "phylum"]:
        for grp, d in pg[unseen].groupby(col):
            if len(d) < 25:
                continue
            r = {"group": f"{col}: {grp}", "n": len(d)}
            r |= {f"KO F1 {METHODS[m][0]}": d[f"ko_f1__{m}"].mean() for m in calls}
            tax.append(r)
    tax = pd.DataFrame(tax)
    if len(tax):
        tax["best"] = tax[[c for c in tax if c.startswith("KO F1")]].idxmax(axis=1).str.removeprefix("KO F1 ")
        tax = tax.sort_values("n", ascending=False)
    tax.to_csv(out / "by_taxon.tsv", sep="\t", index=False)

    u = unseen.to_numpy()
    pk = per_ko(Y[u], {m: s[u] for m, s in scores.items()}, {m: c[u] for m, c in calls.items()}, kos)
    names = pd.read_csv(a.ko_names, sep="\t", header=None, index_col=0).iloc[:, 0]
    pk.insert(1, "name", pk["ko"].map(names))
    pk.to_csv(out / "per_ko.tsv", sep="\t", index=False)
    fig_per_ko(pk, out)
    pk["rf_minus_picrust2"] = pk["auroc__M3_rf"] - pk["auroc__M2_picrust2"]
    top = lambda df: df[["ko", "name", "prevalence", "auroc__M2_picrust2", "auroc__M3_rf", "auroc__M1_knn5"]]  # noqa: E731

    # ---- 3. ordination ----
    ord_methods = ["M2_picrust2", "M1_knn1", "M1_knn5", "M3_rf"]
    rows, stats, coords_all, versus = [], [], [], []
    for label, truth, preds, dist in [
            ("KO presence (Jaccard)", Y[u], {m: calls[m][u] for m in ord_methods}, jaccard_dist),
            ("Module completeness (Bray-Curtis)", Ym[u], {m: mod_pred[m][u] for m in ord_methods}, bray_curtis_dist)]:
        coords, expl, n, st, cr, seg, vs = ordination(truth, preds, dist, label, rng, a.max_ordination)
        versus.append(vs)
        rows.append((label.split(" (")[0], coords, expl, n, st, ord_methods, seg))
        stats.append(st)
        coords_all.append(cr)
    fig_ordination(rows, out)
    stats = pd.concat(stats)
    stats.to_csv(out / "ordination_stats.tsv", sep="\t", index=False)
    versus = pd.concat(versus)
    versus.to_csv(out / "ordination_vs_picrust2.tsv", sep="\t", index=False)
    for df in (stats, versus):  # effect sizes near 0 need 4 decimals in the report tables
        for c in ["PERMANOVA R2", "ANOSIM R"]:
            df[c] = df[c].map(lambda v: f"{v:.4f}")
    co = pd.concat(coords_all)
    co["accession"] = genomes[u][co["genome_index"]]
    co.drop(columns="genome_index").to_csv(out / "ordination_coordinates.tsv", sep="\t", index=False)

    # ---- inspection tables ----
    long = []
    for j, ko in enumerate(kos):
        keep = Y[:, j] | np.any([c[:, j] for c in calls.values()], axis=0)
        if keep.any():
            d = {"accession": genomes[keep], "ko": ko, "truth": Y[keep, j].astype(int),
                 "picrust2_copies": scores["M2_picrust2"][keep, j]}
            d |= {f"{METHODS[m][0]}_prob": scores[m][keep, j].round(3) for m in ours_ko}
            long.append(pd.DataFrame(d))
    pd.concat(long).to_csv(out / "ko_calls.tsv.gz", sep="\t", index=False)
    ml = pd.DataFrame(Ym, index=genomes, columns=mods).stack().rename("truth").to_frame()
    for m, P in mod_pred.items():
        ml[METHODS[m][0]] = pd.DataFrame(P, index=genomes, columns=mods).stack().round(3)
    ml.index.names = ["accession", "module"]
    ml[(ml > 0).any(axis=1)].to_csv(out / "module_calls.tsv.gz", sep="\t")

    ex = pg[unseen].dropna(subset=["nsti"]).copy()
    ex["q"] = pd.qcut(ex["nsti"], 5, labels=False, duplicates="drop")
    ex = ex.groupby("q").apply(lambda d: d.sample(1, random_state=a.seed), include_groups=False)
    ex_rows = []
    for acc, r in ex.droplevel(0).iterrows():
        e = {"accession": acc, "species": r.species, "NSTI": r.nsti, "16S id. to our nearest": r.nn_identity,
             "true KOs": int(r.ko_n_true)}
        for m in ["M2_picrust2", "M1_knn5", "M3_rf"]:
            e[f"{METHODS[m][0]} TP/FP/FN"] = f"{int(r[f'ko_tp__{m}'])}/{int(r[f'ko_fp__{m}'])}/{int(r[f'ko_fn__{m}'])}"
        ex_rows.append(e)

    # ---- report ----
    n_ref = len(genomes)
    sel = lambda df, gs, ms: df[(df.genomes == gs) & df.metric.isin(ms)].drop(columns="genomes")  # noqa: E731
    lines = [
        "# PICRUSt2 against our 16S models", "",
        f"Same genomes, same targets, per genome. {len(genomes)} genomes; **{int(unseen.sum())} are not in "
        f"PICRUSt2's reference** (the fair test set) and {int((~unseen).sum())} are (PICRUSt2 then reads the genome's "
        f"own annotation, a leak). Targets: {len(kos)} KOs (present in 5-95 % of genomes; the Aim 1 set) and "
        f"{len(mods)} KEGG modules. Truth: eggNOG-mapper annotation of each genome.", "",
        "Our predictions are out-of-fold (full-length 16S, random 5-fold split): each genome is predicted by "
        "models trained on the other ~80 % of genomes, never on itself. PICRUSt2 2.6.4 runs with its default "
        f"reference. PICRUSt2 returned no prediction for {len(dropped)} genome(s) (NSTI > 2); they count as "
        "predicting nothing.", "",
        "## 1. Does our approach beat PICRUSt2, and by how much?", "",
        f"**KO sets:** our best model ({best}, mean F1 {b.ours_mean:.3f}) {verdict} PICRUSt2 "
        f"(F1 {b.picrust2_mean:.3f}): difference {b['diff (ours better > 0)']:+.3f} "
        f"(95 % CI {b.ci_lo:+.3f} to {b.ci_hi:+.3f}). Per genome, ours is better on {b.genomes_ours_better:.0%} and "
        f"PICRUSt2 on {b.genomes_picrust2_better:.0%} (ties: the rest): the mean can favour one method while more "
        "genomes favour the other, which is why the Wilcoxon p (ranks) and the CI (mean) can disagree.", "",
        f"**Module completeness:** our best model ({mbest}, mean error {mb.ours_mean:.3f}) {mverdict} PICRUSt2 "
        f"(error {mb.picrust2_mean:.3f}): error lower by {mb['diff (ours better > 0)']:+.3f} "
        f"(95 % CI {mb.ci_lo:+.3f} to {mb.ci_hi:+.3f}).", "",
        "Mean per-genome scores on the unseen genomes. Precision = share of predicted KOs that are real; recall = "
        "share of real KOs that are predicted; `ko_n_predicted` vs `ko_n_true` shows over- or under-prediction.", "",
        md_table(means), "",
        "Paired per genome (our model minus PICRUSt2, signed so that **positive = ours better**; module error "
        "sign-flipped), bootstrap 95 % CI, Wilcoxon signed-rank p, and the share of genomes where ours is better:", "",
        md_table(sel(paired, "not in PICRUSt2 reference", ["ko_f1", "module_mae", "ko_precision", "ko_recall"])), "",
        "For context, the genomes inside PICRUSt2's reference (not a fair test; PICRUSt2 sees their own annotation):",
        "", md_table(sel(paired, "in PICRUSt2 reference (leak)", ["ko_f1"])), "",
        "**Reference size differs.** In this pilot our models learn from ~"
        f"{int(n_ref * 0.8)} genomes per fold; PICRUSt2's reference has 27,870. Accuracy rises with reference "
        "density (Aim 1 figure 1), so this comparison handicaps our models; at full scale ours learn from ~44,000.", "",
        "## 2. Where is either approach better?", "",
        "![novelty](fig1_accuracy_by_novelty.png)", "",
        quote("**Figure 1. Accuracy by how novel each genome is**, unseen genomes only. Left and right: genomes "
              "binned by PICRUSt2's NSTI (tree distance to its nearest reference genome; quintiles). Middle: binned "
              "by 16S identity to our nearest training genome. Each point is the mean per-genome score of one "
              "method in that bin.",
              "**How to read it:** where one line sits above (F1) or below (error) the other, that method is better "
              "for genomes of that novelty. The two x axes measure novelty against different references, so "
              "PICRUSt2's NSTI is the common yardstick for both methods in the left and right panels."), "",
        "By taxon (unseen genomes; groups with >= 25 genomes): mean KO F1 per method.", "",
        md_table(tax) if len(tax) else "_no group with >= 25 genomes_", "",
        "![per KO](fig2_per_ko_auroc.png)", "",
        quote("**Figure 2. Per-KO ranking quality**, unseen genomes. Each dot is one KO; AUROC measures how well a "
              "method's score (PICRUSt2: predicted copy number; ours: probability) ranks genomes that carry the KO "
              "above genomes that do not (0.5 = chance, 1 = perfect).",
              "Right: per-KO F1 of the presence calls (PICRUSt2: copy number >= 1; 5-NN: probability >= 0.5).",
              "**How to read it:** dots above the diagonal are KOs our model predicts better; below, PICRUSt2 "
              "does. AUROC rewards a graded score: our probabilities are graded, PICRUSt2's copy numbers are coarse "
              "(mostly 0, 1, 2), which limits its AUROC. The F1 panel compares the yes/no calls that users read.",
              f"**Here:** median per-KO AUROC {pk['auroc__M3_rf'].median():.3f} (random forest) vs "
              f"{pk['auroc__M2_picrust2'].median():.3f} (PICRUSt2); median per-KO F1 "
              f"{pk['f1__M1_knn5'].median():.3f} (5-NN) vs {pk['f1__M2_picrust2'].median():.3f} (PICRUSt2)."), "",
        "KOs where the random forest beats PICRUSt2 most (AUROC):", "",
        md_table(top(pk.nlargest(10, "rf_minus_picrust2"))), "",
        "KOs where PICRUSt2 beats the random forest most:", "",
        md_table(top(pk.nsmallest(10, "rf_minus_picrust2"))), "",
        "## 3. Ordination: predicted against true functional profiles", "",
        "![ordination](fig3_ordination.png)", "",
        quote("**Figure 3. Do the predicted profiles reproduce the functional landscape of the genomes?** For each "
              "row, the true profiles of the unseen genomes (grey) and every method's predicted profiles for the "
              "same genomes are placed in one principal-coordinates analysis (PCoA), so all panels of a row share "
              "axes. Top: KO presence/absence, Jaccard distance. Bottom: KEGG module completeness, Bray-Curtis "
              "distance. Thin lines join a genome's true profile to its predicted profile (random subset of genomes)."
              + (f" Ordination uses a random {a.max_ordination} genomes." if unseen.sum() > a.max_ordination else ""),
              "**How to read it:** a good method's coloured cloud overlays the grey cloud and its lines are short. "
              "Each panel's box gives two kinds of numbers. *Agreement, genome by genome:* Mantel r = correlation "
              "between true and predicted genome-to-genome distances (1 = the method keeps which genomes are "
              "functionally alike); Procrustes m² = mismatch after best rotation and scaling of the two 10-axis "
              "ordinations (0 = identical). *Group tests, predicted cloud against true cloud:* PERMANOVA R² = share "
              "of all variation explained by 'predicted vs true' (0 = clouds share a centre); ANOSIM R = how much "
              "farther apart profiles are across the two groups than within (0 = no separation, 1 = complete); "
              "spread = mean distance to the group centroid relative to truth (PERMDISP; below 1 = predictions are "
              "more alike than real genomes, i.e. pulled towards typical profiles).",
              "Tests use 999 permutations that swap the two labels within each genome, because every genome appears "
              "in both groups (paired design). With about a thousand genomes, very small differences reach "
              "p = 0.001: read the effect sizes (R², R, spread), not only p. Under these within-genome swaps "
              "the null ANOSIM R sits slightly below 0, so an R of about 0 can still give p = 0.001."), "",
        "Agreement with the truth, genome by genome:", "",
        md_table(stats[["space", "method", "genomes", "mantel_r (between-genome distances vs truth)",
                        "procrustes_m2 (0 = same configuration)", "mean distance to own truth"]]), "",
        "Group tests, each method's predicted profiles against the true profiles:", "",
        md_table(stats[["space", "method", "PERMANOVA R2", "PERMANOVA p", "ANOSIM R", "ANOSIM p", "spread truth",
                        "spread method", "PERMDISP p"]]), "",
        "Group tests, our predicted profiles against PICRUSt2's (same genomes):", "",
        md_table(versus), "",
        "## 4. How the methods differ", "",
        md_table(pd.DataFrame([(k, o.format(n_ref=n_ref), p) for k, o, p in DIFFERENCES],
                              columns=["", "Our approach", "PICRUSt2 2.6.4 (default)"])), "",
        "PICRUSt2 facts: Douglas et al. 2020 (Nat Biotechnol 38:685); PICRUSt2-SC reference database "
        "(Bioinformatics 2025, btaf269); genome counts read from the installed reference.", "",
        "## 5. Inspect the outputs", "",
        "All in the run's results folder (`results/<profile>/`):", "",
        "- `picrust2_comparison/per_genome.tsv`: one row per genome: taxonomy, in-reference flag, NSTI, our nearest "
        "16S identity, and per method KO F1, precision, recall, TP/FP/FN, module error.",
        "- `picrust2_comparison/ko_calls.tsv.gz`: genome x KO, every pair that is true or predicted by any method: "
        "truth, PICRUSt2 copy number, our probabilities.",
        "- `picrust2_comparison/module_calls.tsv.gz`: genome x module completeness: truth and every method.",
        "- `picrust2_comparison/per_ko.tsv`: per KO, with name: prevalence, AUROC and F1 per method.",
        "- `picrust2_comparison/paired.tsv`, `by_novelty.tsv`, `by_taxon.tsv`, `ordination_stats.tsv`, "
        "`ordination_vs_picrust2.tsv`, "
        "`ordination_coordinates.tsv`: the numbers behind every table and figure.",
        "- `aim1/<arm>_<level>/predictions/<model>_<target>.parquet`: all out-of-fold predictions of our models, "
        "every arm and holdout level.",
        "- `picrust2_out/`: PICRUSt2's own outputs (KO copy numbers, NSTI, reference genome list).",
        "- `targets/`: the truth (KO presence, module completeness, COG profiles).", "",
        "Example genomes, one per NSTI quintile (unseen; TP = true KOs predicted, FP = predicted but absent, "
        "FN = present but missed):", "",
        md_table(pd.DataFrame(ex_rows)), "",
    ]
    (out / "picrust2_comparison.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

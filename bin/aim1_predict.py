#!/usr/bin/env python3
"""Aim 1: predict functional targets from 16S allele features, for one input arm and one holdout level.

Usage: aim1_predict.py --features features/ --targets targets/ --splits splits.tsv --arm full --level family
                       [--ko-min-prev 0.05] [--ko-max-prev 0.95] [--threads 8] [--seed 42] --outdir out

Models (DESIGN §4.1): M0 null (training prevalence/mean); M1 k-NN on 16S identity, k=1 and
distance-weighted k=5; M3 multi-output random forest (all target types); M4 LightGBM per target
(modules and COG only). M2 (PICRUSt2) runs in its own step.

Targets are filtered once on all analysed genomes: KO prevalence in [min, max]; modules with
completeness > 0 in >= 5 % of genomes and variance > 1e-3; all COG categories.

Outputs (outdir):
  genome_metrics.tsv  accession, fold, model, target_type, metric, value   (per test genome)
  target_metrics.tsv  target, target_type, model, metric, value            (pooled over folds)
  nearest.tsv         accession, fold, nn_identity                         (16S identity to closest training genome)
  predictions/<model>_<target_type>.parquet   out-of-fold predictions, genome x target (KO: probability; modules: completeness; COG: fraction)
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import average_precision_score, roc_auc_score

BASES = np.array(list("ACGT"))


# ---------- data ----------

def arm_columns(features_dir, arm):
    """Alignment columns covered by an arm (full = all kept columns)."""
    cols = pd.read_csv(Path(features_dir) / "alignment_columns.tsv", sep="\t")["column"]
    if arm == "full":
        return cols.astype(str).tolist()
    meta = pd.read_csv(Path(features_dir) / "feature_columns.tsv", sep="\t")
    inside = meta.loc[meta[f"in_{arm}"] == 1, "column"]
    return cols[cols.between(inside.min(), inside.max())].astype(str).tolist()


def load(features_dir, targets_dir, arm, ko_min, ko_max):
    X = pd.read_parquet(Path(features_dir) / f"features_{arm}.parquet")
    chars = pd.read_parquet(Path(features_dir) / "alignment_kept.parquet", columns=arm_columns(features_dir, arm))
    t = Path(targets_dir)
    ko = pd.read_parquet(t / "ko_presence.parquet")
    mod = pd.read_parquet(t / "module_completeness.parquet")
    cog = pd.read_parquet(t / "cog_profile.parquet")
    genomes = X.index.intersection(ko.index).sort_values()
    ko, mod, cog = ko.loc[genomes], mod.loc[genomes], cog.loc[genomes]
    prev = ko.mean()
    ko = ko.loc[:, prev.between(ko_min, ko_max)]
    mod = mod.loc[:, ((mod > 0).mean() >= 0.05) & (mod.var() > 1e-3)]
    targets = {"ko": ko.astype("float32"), "module": mod.astype("float32"), "cog": cog.astype("float32")}
    return X.loc[genomes], chars.loc[genomes], targets


def base_onehot(chars):
    """genome x (column, base) 0/1 matrix; gaps and ambiguity codes are all-zero."""
    a = chars.to_numpy().astype("U1")
    return (a[:, :, None] == BASES[None, None, :]).reshape(len(a), -1).astype("float32")


def identity(onehot_test, onehot_train):
    """Matching bases / bases in the longer of the two sequences (over the arm's columns)."""
    m = onehot_test @ onehot_train.T
    n_test, n_train = onehot_test.sum(1), onehot_train.sum(1)
    return m / np.maximum(n_test[:, None], n_train[None, :])


# ---------- models ----------

def knn_predict(sim, Y_train, k):
    if k == 1:
        return Y_train[np.argmax(sim, axis=1)]
    idx = np.argsort(-sim, axis=1)[:, :k]
    w = 1.0 / (1.0 - np.take_along_axis(sim, idx, axis=1) + 1e-3)
    w /= w.sum(1, keepdims=True)
    return np.einsum("nk,nkt->nt", w, Y_train[idx])


def rf_predict(X_train, Y_train, X_test, threads, seed):
    # Multi-output regression on 0/1 or fractions: one forest for all targets of a type; leaf
    # means are probabilities for KOs. Fixed hyperparameters in the pilot (DECISIONS D048).
    rf = RandomForestRegressor(n_estimators=100, max_features="sqrt", min_samples_leaf=3,
                               n_jobs=threads, random_state=seed)
    rf.fit(X_train, Y_train)
    return rf.predict(X_test).reshape(len(X_test), -1)


def lgbm_predict(X_train, Y_train, X_test, threads, seed):
    def one(y):
        if np.ptp(y) == 0:
            return np.full(len(X_test), y[0])
        m = LGBMRegressor(n_estimators=200, learning_rate=0.05, num_leaves=15, min_child_samples=10,
                          colsample_bytree=0.5, subsample=0.8, subsample_freq=1, n_jobs=1,
                          random_state=seed, verbose=-1)
        return m.fit(X_train, y).predict(X_test)
    cols = Parallel(n_jobs=threads)(delayed(one)(Y_train[:, j]) for j in range(Y_train.shape[1]))
    return np.column_stack(cols)


# ---------- metrics ----------

def genome_metrics(ttype, Y, P):
    """Per test genome: KO set F1/Jaccard at 0.5; module MAE/Spearman; COG Spearman/Bray-Curtis."""
    out = {}
    if ttype == "ko":
        pred, true = P >= 0.5, Y >= 0.5
        tp = (pred & true).sum(1)
        fp = (pred & ~true).sum(1)
        fn = (~pred & true).sum(1)
        out["f1"] = 2 * tp / np.maximum(2 * tp + fp + fn, 1)
        out["jaccard"] = tp / np.maximum(tp + fp + fn, 1)
    else:
        out["mae"] = np.abs(P - Y).mean(1)
        out["spearman"] = np.array([spearmanr(p, y).statistic if np.ptp(y) and np.ptp(p) else np.nan
                                    for p, y in zip(P, Y)])
        if ttype == "cog":
            out["bray_curtis"] = np.abs(P - Y).sum(1) / np.maximum((P + Y).sum(1), 1e-12)
    return out


def target_metrics(ttype, Y, P):
    """Per target, pooled over all test folds."""
    rows = {}
    for j in range(Y.shape[1]):
        y, p = Y[:, j], P[:, j]
        if ttype == "ko":
            if 0 < y.sum() < len(y):
                rows[j] = {"auroc": roc_auc_score(y, p), "auprc": average_precision_score(y, p),
                           "prevalence": y.mean()}
        else:
            rows[j] = {"spearman": spearmanr(p, y).statistic if np.ptp(p) and np.ptp(y) else np.nan,
                       "mae": np.abs(p - y).mean()}
    return rows


# ---------- main ----------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--features", required=True)
    ap.add_argument("--targets", required=True)
    ap.add_argument("--splits", required=True)
    ap.add_argument("--arm", required=True)
    ap.add_argument("--level", required=True)
    ap.add_argument("--ko-min-prev", type=float, default=0.05)
    ap.add_argument("--ko-max-prev", type=float, default=0.95)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    X, chars, targets = load(a.features, a.targets, a.arm, a.ko_min_prev, a.ko_max_prev)
    folds = pd.read_csv(a.splits, sep="\t").query("level == @a.level").set_index("accession")["fold"]
    folds = folds.reindex(X.index).dropna().astype(int)
    X, chars = X.loc[folds.index], chars.loc[folds.index]
    Xv, oh = X.to_numpy(), base_onehot(chars)
    n = len(X)

    preds = {}  # (model, ttype) -> n x targets
    nearest = []
    for f in sorted(folds.unique()):
        te = np.flatnonzero(folds.to_numpy() == f)
        tr = np.flatnonzero(folds.to_numpy() != f)
        sim = identity(oh[te], oh[tr])
        nearest.append(pd.DataFrame({"accession": X.index[te], "fold": f, "nn_identity": sim.max(1)}))
        for ttype, T in targets.items():
            Y = T.to_numpy()
            fold_preds = {
                "M0_null": np.tile(Y[tr].mean(0), (len(te), 1)),
                "M1_knn1": knn_predict(sim, Y[tr], 1),
                "M1_knn5": knn_predict(sim, Y[tr], 5),
                "M3_rf": rf_predict(Xv[tr], Y[tr], Xv[te], a.threads, a.seed),
            }
            if ttype != "ko":
                fold_preds["M4_lgbm"] = lgbm_predict(Xv[tr], Y[tr], Xv[te], a.threads, a.seed)
            for model, P in fold_preds.items():
                preds.setdefault((model, ttype), np.full((n, Y.shape[1]), np.nan, dtype="float32"))[te] = P

    g_rows, t_rows = [], []
    for (model, ttype), P in preds.items():
        T = targets[ttype]
        Y = T.to_numpy()
        for metric, vals in genome_metrics(ttype, Y, P).items():
            g_rows.append(pd.DataFrame({"accession": X.index, "fold": folds.to_numpy(), "model": model,
                                        "target_type": ttype, "metric": metric, "value": vals}))
        for j, ms in target_metrics(ttype, Y, P).items():
            t_rows += [{"target": T.columns[j], "target_type": ttype, "model": model, "metric": k, "value": v}
                       for k, v in ms.items()]

    # out-of-fold predictions, for inspection and the PICRUSt2 comparison: one genome x target table per model
    (out / "predictions").mkdir(exist_ok=True)
    for (model, ttype), P in preds.items():
        pd.DataFrame(P, index=X.index, columns=targets[ttype].columns).to_parquet(
            out / "predictions" / f"{model}_{ttype}.parquet", compression="zstd")

    meta = {"arm": a.arm, "level": a.level}
    pd.concat(g_rows).assign(**meta).to_csv(out / "genome_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(t_rows).assign(**meta).to_csv(out / "target_metrics.tsv", sep="\t", index=False)
    pd.concat(nearest).assign(**meta).to_csv(out / "nearest.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()

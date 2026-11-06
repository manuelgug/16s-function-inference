#!/usr/bin/env python3
"""Exploration harness: candidate ideas for community function profiles from ASVs, scored against PICRUSt2.

Usage (inside the analysis image, from the repo root):
  python exploration/explore.py --data results/bench/bench_data --tasks results/bench/tasks/*_0_* \
      --replicates 5 --outdir results/exploration

Leakage rules (DESIGN §B, D095):
- Only dev-partition genomes are queried (the task dirs come from a BENCH_PHASE=dev run); confirm genomes are
  never touched.
- References = the task's training genomes (held-out clade removed), exactly as for PICRUSt2 in the same task.
- Every learned component is cross-fitted by family inside the task: fitted on the dev queries of one half of
  the families, applied to the other half, and vice versa. No genome is scored by a model that saw its family.
- Communities are simulated only from the task's held-out genomes (bin/community.py).
Every variant and its result is written to results.tsv; exploration/report.py turns it into the report.
"""
import argparse
import sys
import tempfile
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import s16fun as sf  # noqa: E402
from benchlib import WeightedAUROC, ko_scores, load_data, load_pred, write_fasta  # noqa: E402
from community import HABITATS, profile, score, simulate  # noqa: E402

P2 = ["mp", "emp_prob", "pic", "subtree_average"]
LGB = dict(learning_rate=0.1, num_leaves=63, min_data_in_leaf=200, verbose=-1)


def crossfit(X_rows, y_rows, groups, params, rounds, threads, rng, max_rows=4_000_000):
    """Fit on one family half, predict the other (and vice versa). X_rows: callable(rows) -> feature matrix."""
    preds = {}
    for g in (0, 1):
        tr, te = np.flatnonzero(groups != g), np.flatnonzero(groups == g)
        X, y = X_rows(tr), y_rows(tr)
        sel = rng.choice(len(y), min(len(y), max_rows), replace=False)
        m = lgb.train({**params, "num_threads": threads, "seed": int(rng.integers(1 << 31))},
                      lgb.Dataset(X[sel], y[sel]), rounds)
        preds[g] = (te, m)
    return preds


def ratio_evidence(ref, nbr, w):
    """Neighbour-weighted KO copies per 16S copy: mean over the 25 nearest genomes of C_jk / s_j (all of them),
    and over those carrying the KO."""
    s = ref.genomes.ssu_copies.to_numpy(float)
    R = sparse.diags(1.0 / s) @ ref.copies.astype(np.float32)
    R = sparse.csr_matrix(R)
    all_ = sf.knn_average(R, nbr, w, sf.K)
    carry = sf.knn_average(ref.Y, nbr, w, sf.K)
    return all_, np.where(carry > 0, all_ / np.maximum(carry, 1e-12), 0.0)


def run_task(t, data, threads, replicates, rng):
    genomes, kos, copies, seqs = data
    info = pd.read_csv(t / "task.tsv", sep="\t")
    level, arm = info.level.iloc[0], info.arm.iloc[0]
    q = info.accession.tolist()
    full = sf.Reference(genomes, None, kos, copies)
    train = (t / "train.txt").read_text().split()
    with tempfile.TemporaryDirectory() as tmp:
        write_fasta(Path(tmp) / "ref.fna", seqs, train)
        ref = full.subset(train, Path(tmp) / "ref.fna")
        hits = sf.search(t / "query.fna", ref.ssu_fna, threads, tmp)
    ev, nov = sf.evidence(hits, q, ref)
    nbr, ids, w = sf.neighbours(hits, q, ref)
    row = [full.row[x] for x in q]
    Y = full.Y[row].toarray() > 0
    T = np.minimum(full.copies[row].toarray(), 10).astype(np.float32)
    T16 = genomes.loc[q, "ssu_copies"].to_numpy(float)
    fam = genomes.loc[q, "family"].to_numpy()
    uf = np.unique(fam)
    half = dict(zip(uf, rng.permutation(len(uf)) % 2))
    groups = np.array([half[f] for f in fam])
    nk = len(kos)

    # ---- learned pieces, cross-fitted by family ----
    P = np.zeros((len(q), nk), np.float32)
    for g, (te, m) in crossfit(lambda r: sf.feature_rows(ev, r), lambda r: Y[r].ravel(), groups,
                               {**LGB, "objective": "binary"}, 300, threads, rng).items():
        P[te] = sf.probabilities({k: v[te] for k, v in ev.items()}, m)
    ssu_l = np.ones(len(q))
    for g, (te, m) in crossfit(lambda r: ev["_ssu_X"][r], lambda r: np.log(T16[r]), groups,
                               dict(objective="l1", learning_rate=0.05, num_leaves=15, min_data_in_leaf=20, verbose=-1),
                               400, threads, rng).items():
        ssu_l[te] = np.maximum(np.exp(m.predict(ev["_ssu_X"][te])), 1)
    cond = ev["_cond_copies"]
    knn_copies = sf.knn_average(ref.copies.astype(np.float32).tocsr(), nbr, w, 5)
    ev_cp = dict(ev, knn5_copies=knn_copies, cond_copies=cond)

    def regress(evd, feats, target):
        """Cross-fitted Tweedie regression of a per genome x KO target on evidence features."""
        rows_X = lambda r: np.stack([np.asarray(evd[f][r], np.float32).ravel() for f in feats], 1)  # noqa: E731
        out = np.zeros((len(q), nk), np.float32)
        for g, (te, m) in crossfit(rows_X, lambda r: target[r].ravel(), groups,
                                   {**LGB, "objective": "tweedie", "tweedie_variance_power": 1.2}, 300, threads, rng).items():
            out[te] = np.maximum(m.predict(rows_X(te)), 0).reshape(len(te), nk)
        return out

    Creg = regress(ev_cp, sf.FEATURES + ["knn5_copies", "cond_copies"], T)
    ratio_all, ratio_carry = ratio_evidence(ref, nbr, w)
    Rreg = regress(dict(ev, ratio_all=ratio_all, ratio_carry=ratio_carry), sf.FEATURES + ["ratio_all", "ratio_carry"],
                   T / T16[:, None])
    median16 = nov.ssu_copies.to_numpy()

    # ---- PICRUSt2 ----
    p2 = {m: load_pred(t / f"pred_picrust2_{m}.npz") for m in P2}
    for m in P2:
        assert p2[m]["queries"] == q, "PICRUSt2 queries differ"
    p2_16 = pd.read_csv(t / "nsti.tsv", sep="\t", index_col=0).reindex(q).ssu_copies.fillna(1).to_numpy(float)

    # ---- variants: (expected KO copies per genome, 16S copies) or ratio (KO copies per 16S copy) ----
    V = {f"PICRUSt2 {m}": dict(copies=p2[m]["copies"], ssu=p2_16, score=p2[m]["score"], call=p2[m]["call"]) for m in P2}
    knn16 = np.maximum(sf.knn_average(sparse.csr_matrix(ref.genomes[["ssu_copies"]].to_numpy(np.float32)), nbr, w, 5)[:, 0], 1)
    V["E1 5-NN copies"] = dict(copies=knn_copies, ssu=knn16, score=ev["knn5"], call=ev["knn5"] >= 0.5)
    V["E2 s16fun p x copies, median 16S"] = dict(copies=P * cond, ssu=median16, score=P, call=sf.expected_f1_calls(P))
    V["E3 s16fun p x copies, learned 16S"] = dict(copies=P * cond, ssu=ssu_l, score=P, call=V["E2 s16fun p x copies, median 16S"]["call"])
    V["E3h s16fun hard calls, learned 16S"] = dict(copies=V["E3 s16fun p x copies, learned 16S"]["call"] * cond, ssu=ssu_l)
    V["E4 copy regression (Tweedie), learned 16S"] = dict(copies=Creg, ssu=ssu_l, score=Creg)
    V["E5 ratio, 25-NN KO copies per 16S copy"] = dict(ratio=ratio_all)
    V["E5b ratio, p x ratio when present"] = dict(ratio=P * ratio_carry)
    V["E6 ratio regression (Tweedie)"] = dict(ratio=Rreg)
    V["D1 s16fun with TRUE 16S copies"] = dict(copies=P * cond, ssu=T16)
    V["D2 s16fun KO + PICRUSt2 16S copies"] = dict(copies=P * cond, ssu=p2_16)
    best_p2 = "pic"
    V["D3 mean of s16fun and PICRUSt2 pic"] = dict(copies=(P * cond + p2[best_p2]["copies"]) / 2, ssu=np.sqrt(ssu_l * p2_16))

    rows = []
    tag = dict(level=level, arm=arm, genomes=len(q))
    # per-genome metrics
    for name, v in V.items():
        if "call" in v:
            rows.append(dict(**tag, variant=name, metric="ko_f1", value=float(ko_scores(Y, v["call"])["ko_f1"].mean())))
        if "copies" in v:
            rows.append(dict(**tag, variant=name, metric="copy_mae",
                             value=float(np.abs(np.minimum(v["copies"], 10) - T).mean())))
        if "ssu" in v:
            rows.append(dict(**tag, variant=name, metric="ssu_log_error",
                             value=float(np.abs(np.log(np.maximum(v["ssu"], 1) / T16)).mean())))
        if "score" in v:
            keep = (Y.sum(0) >= 5) & ((~Y).sum(0) >= 5)
            a = WeightedAUROC(v["score"][:, keep], Y[:, keep])(np.ones(len(q)))
            rows.append(dict(**tag, variant=name, metric="ko_auroc_mean", value=float(np.nanmean(a))))
    # community metrics
    for habitat in HABITATS:
        for rep in range(replicates):
            cells, group, _ = simulate(len(q), habitat, rng)
            reads = cells * T16
            truth = profile(cells, T)
            for name, v in V.items():
                X = profile(reads, v["ratio"]) if "ratio" in v else profile(reads / np.maximum(v["ssu"], 1), v["copies"])
                for k, val in score(X, truth, group).items():
                    rows.append(dict(**tag, variant=name, metric=f"community_{k}", habitat=habitat, replicate=rep, value=val))
    print(f"{t.name}: {len(q)} genomes", flush=True)
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True)
    ap.add_argument("--tasks", nargs="+", required=True)
    ap.add_argument("--replicates", type=int, default=5)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    data = load_data(a.data)
    rng = np.random.default_rng(a.seed)
    rows = []
    for t in sorted(Path(x) for x in a.tasks):
        info = pd.read_csv(t / "task.tsv", sep="\t")
        assert (info.phase == "dev").all(), f"{t}: exploration uses dev genomes only"
        rows += run_task(t, data, a.threads, a.replicates, rng)
        pd.DataFrame(rows).to_csv(out / "results.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()

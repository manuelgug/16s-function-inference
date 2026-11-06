#!/usr/bin/env python3
"""Our method and the baseline methods on one benchmark task (output of bench_prepare.py).

Usage: bench_predict.py --data bench_data/ --task task/ [--stackers stackers/] [--sample-genomes 150]
                        [--threads 16] [--seed 42] --outdir pred/

Without --stackers (pass 1): baselines that need no training beyond the reference, plus a sample of
stacker training rows (dev phase only):
  null         KO prevalence in the training genomes (call: >= 0.5)
  nn1          nearest 16S hit's KOs, no identity cutoff (k = 1 nearest neighbour)
  knn5         identity-weighted average of the 5 nearest genomes (call: >= 0.5)
  nearest97    nearest hit's KOs if its identity >= 97 %, else nothing (Piphillin / Tax4Fun2 family)
  genuscore97  KOs in >= 50 % of the nearest hit's genus, if identity >= 97 %, else nothing
               (MicFunPred / PanFP family: genus core genome)
  oracle5      NOT a method: majority vote of the 5 training genomes whose KO sets are most similar to the
               query's true KO set (uses the truth). The ceiling for any method that predicts from relatives.
With --stackers (pass 2), our method and its ablations (and novelty.tsv with the learned 16S copy numbers):
  s16fun       stacked probabilities, expected-F1 calls (the method as shipped, bin/s16fun.py)
  s16fun_t05   same probabilities, calls at >= 0.5 (ablation: decision rule)
  knn5_ef1     5-NN average with expected-F1 calls (ablation: no stacker)
Dev queries use the cross-fitted stacker of their group (never trained on their own family);
confirm queries use the stacker trained on all dev rows.

Outputs: pred_<method>.npz (benchlib.save_pred), novelty.tsv, [stacker_rows.npz, ssu_rows.npz: X, y, accession].
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import s16fun as sf
from benchlib import load_data, save_pred


def oracle5(Yq, Yt, chunk=1000):
    """Mean KO presence of the 5 training genomes with the highest Jaccard similarity to each query's truth."""
    nt = np.asarray(Yt.sum(1)).ravel()
    out = np.zeros((Yq.shape[0], Yt.shape[1]), np.float32)
    for s in range(0, Yq.shape[0], chunk):
        q = Yq[s:s + chunk]
        inter = (q @ Yt.T).toarray()
        jac = inter / np.maximum(np.asarray(q.sum(1)) + nt[None, :] - inter, 1)
        top = np.argpartition(-jac, 4, axis=1)[:, :5]
        out[s:s + chunk] = Yt[top.ravel()].toarray().reshape(len(top), 5, -1).mean(1)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True)
    ap.add_argument("--task", required=True)
    ap.add_argument("--stackers")
    ap.add_argument("--sample-genomes", type=int, default=150)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    out, task = Path(a.outdir), Path(a.task)
    out.mkdir(parents=True, exist_ok=True)

    genomes, kos, copies, _ = load_data(a.data)
    full = sf.Reference(genomes, None, kos, copies)
    train = (task / "train.txt").read_text().split()
    ref = full.subset(train, task / "ref.fna")
    info = pd.read_csv(task / "task.tsv", sep="\t")
    queries = info.accession.tolist()
    hits = sf.search(task / "query.fna", ref.ssu_fna, a.threads, out)
    (out / "hits.tsv").unlink()
    ev, nov = sf.evidence(hits, queries, ref)
    nov.to_csv(out / "novelty.tsv", sep="\t")
    top = nov.top_id.to_numpy()[:, None]

    if not a.stackers:
        near = top >= 0.97
        base = {
            "null": (ev["prevalence"], ev["prevalence"] >= 0.5),
            "nn1": (ev["knn1"], ev["knn1"] >= 0.5),
            "knn5": (ev["knn5"], ev["knn5"] >= 0.5),
            "nearest97": (ev["knn1"] * near, (ev["knn1"] >= 0.5) & near),
            "genuscore97": (ev["prior_genus"] * near, (ev["prior_genus"] >= 0.5) & near),
        }
        for m, (score, call) in base.items():
            save_pred(out / f"pred_{m}.npz", queries, score, call)
        O = oracle5(full.Y[[full.row[q] for q in queries]], ref.Y)
        save_pred(out / "pred_oracle5.npz", queries, O, O >= 0.5)
        if info.phase.iloc[0] == "dev" and queries:
            rng = np.random.default_rng(a.seed + int(info.fold.iloc[0]))
            rows = np.sort(rng.choice(len(queries), min(a.sample_genomes, len(queries)), replace=False))
            Y = full.Y[[full.row[q] for q in np.array(queries)[rows]]].toarray() > 0
            np.savez_compressed(out / "stacker_rows.npz", X=sf.feature_rows(ev, rows).astype(np.float16),
                                y=Y.ravel(), accession=np.repeat(np.array(queries)[rows], len(kos)))
            np.savez_compressed(out / "ssu_rows.npz", X=ev["_ssu_X"],
                                y=np.log(genomes.loc[queries, "ssu_copies"].to_numpy(float)), accession=np.array(queries))
        return

    st = Path(a.stackers)
    P = np.zeros((len(queries), len(kos)), np.float32)
    ssu = np.ones(len(queries))
    if info.phase.iloc[0] == "dev":
        group = pd.read_csv(st / "crossfit.tsv", sep="\t", index_col=0)["group"]
        for gid in sorted(set(group.reindex(queries))):
            rows = np.flatnonzero(group.reindex(queries).to_numpy() == gid)
            sub = {k: v[rows] for k, v in ev.items()}
            P[rows] = sf.probabilities(sub, sf.load_stacker(st / f"stacker_g{gid}.txt"))
            ssu[rows] = sf.ssu_copies(sub, nov.iloc[rows], sf.load_stacker(st / f"ssu_g{gid}.txt"))
    else:
        P = sf.probabilities(ev, sf.load_stacker(st / "stacker_all.txt"))
        ssu = sf.ssu_copies(ev, nov, sf.load_stacker(st / "ssu_all.txt"))
    nov.assign(ssu_copies=ssu).to_csv(out / "novelty.tsv", sep="\t")
    save_pred(out / "pred_s16fun.npz", queries, P, sf.expected_f1_calls(P), copies=P * ev["_cond_copies"])
    save_pred(out / "pred_s16fun_t05.npz", queries, P, P >= 0.5)
    save_pred(out / "pred_knn5_ef1.npz", queries, ev["knn5"], sf.expected_f1_calls(ev["knn5"]))


if __name__ == "__main__":
    main()

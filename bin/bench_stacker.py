#!/usr/bin/env python3
"""Fit the s16fun stacker on dev-phase evidence rows, cross-fitted by family.

Usage: bench_stacker.py --rows stacker_rows.npz [...] --ssu-rows ssu_rows.npz [...] --genomes genomes.tsv
                        --splits splits.tsv [--groups 5]
                        [--max-rows 20000000] [--rounds 300] [--threads 16] [--seed 42] --outdir stackers/

Rows come from bench_predict.py pass 1 on dev tasks (every level and arm): evidence features of
sampled held-out genomes x every KO, with the true presence. Dev families are split into --groups
groups; stacker_g<k> is trained without group k's genomes and scores only group k's genomes in the dev
evaluation, so no dev score comes from a model that saw that genome's family. stacker_all is trained
on every dev row; it is the frozen model for the confirm run and the model the tool ships.

The 16S copy-number model (ssu_g<k>, ssu_all) is fitted the same way on ssu_rows (one row per dev genome:
copy numbers and identities of its nearest training genomes; target log copies, L1 loss).

Outputs: crossfit.tsv (accession, group), stacker_g<k>.txt, stacker_all.txt, ssu_g<k>.txt, ssu_all.txt
(LightGBM text models), stacker_importance.tsv.
"""
import argparse
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd

from s16fun import FEATURES

PARAMS = dict(objective="binary", learning_rate=0.1, num_leaves=63, min_data_in_leaf=200, verbose=-1)
SSU_PARAMS = dict(objective="l1", learning_rate=0.05, num_leaves=15, min_data_in_leaf=20, verbose=-1)


def fit(X, y, rounds, threads, rng, max_rows):
    sel = rng.choice(len(y), min(len(y), max_rows), replace=False) if len(y) > max_rows else slice(None)
    return lgb.train({**PARAMS, "num_threads": threads, "seed": int(rng.integers(1 << 31))},
                     lgb.Dataset(X[sel], y[sel], feature_name=FEATURES), rounds)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rows", nargs="+", required=True)
    ap.add_argument("--ssu-rows", nargs="+", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--splits", required=True)
    ap.add_argument("--groups", type=int, default=5)
    ap.add_argument("--max-rows", type=int, default=20_000_000)
    ap.add_argument("--rounds", type=int, default=300)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    parts = [np.load(p) for p in sorted(a.rows)]
    X = np.concatenate([p["X"] for p in parts]).astype(np.float32)
    y = np.concatenate([p["y"] for p in parts])
    acc = np.concatenate([p["accession"] for p in parts])

    g = pd.read_csv(a.genomes, sep="\t", index_col=0)
    s = pd.read_csv(a.splits, sep="\t").drop_duplicates("accession").set_index("accession")
    dev = s.index[s.partition == "dev"].sort_values()
    fams = np.array(sorted(set(g.loc[dev, "family"])), dtype=object)
    fam_group = dict(zip(fams[rng.permutation(len(fams))], np.arange(len(fams)) % a.groups))
    cf = pd.DataFrame({"group": [fam_group[f] for f in g.loc[dev, "family"]]}, index=dev)
    cf.to_csv(out / "crossfit.tsv", sep="\t")
    row_group = cf.loc[acc, "group"].to_numpy()

    for k in range(a.groups):
        fit(X[row_group != k], y[row_group != k], a.rounds, a.threads, rng, a.max_rows).save_model(
            str(out / f"stacker_g{k}.txt"))
    m = fit(X, y, a.rounds, a.threads, rng, a.max_rows)
    m.save_model(str(out / "stacker_all.txt"))

    ssu = [np.load(p) for p in sorted(a.ssu_rows)]
    sX = np.concatenate([p["X"] for p in ssu])
    sy = np.concatenate([p["y"] for p in ssu])
    sg = cf.loc[np.concatenate([p["accession"] for p in ssu]), "group"].to_numpy()
    for k in list(range(a.groups)) + ["all"]:
        sel = sg != k if k != "all" else slice(None)
        lgb.train({**SSU_PARAMS, "num_threads": a.threads, "seed": int(rng.integers(1 << 31))},
                  lgb.Dataset(sX[sel], sy[sel]), 400).save_model(str(out / f"ssu_{'g' + str(k) if k != 'all' else 'all'}.txt"))
    imp = pd.Series(m.feature_importance("gain"), index=FEATURES)
    (imp / imp.sum()).rename("gain_share").to_csv(out / "stacker_importance.tsv", sep="\t")
    print(f"{len(y)} rows from {len(np.unique(acc))} genomes; {a.groups} cross-fit groups")


if __name__ == "__main__":
    main()

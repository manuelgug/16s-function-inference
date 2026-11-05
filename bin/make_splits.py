#!/usr/bin/env python3
"""Cross-validation folds at nested phylogenetic holdout levels (DESIGN §4.2).

Usage: make_splits.py --genomes genomes.tsv --usable genome_qc.tsv --levels random,genus,family,order,class
                      --folds 5 --seed 42 --out splits.tsv

level "random": genomes assigned to folds at random (species-level holdout, since each genome is
its own GTDB species). Other levels: whole taxa at that rank are held out together (GroupKFold),
so no test genome shares that taxon with any training genome.
Output: long table accession, level, fold. Only usable genomes (genome_qc.tsv usable == True) are split.
"""
import argparse

import pandas as pd
from sklearn.model_selection import GroupKFold, KFold

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def make_splits(df, levels, folds, seed):
    out = []
    for level in levels:
        if level == "random":
            cv = KFold(n_splits=folds, shuffle=True, random_state=seed).split(df)
        elif level in RANKS:
            cv = GroupKFold(n_splits=folds, shuffle=True, random_state=seed).split(df, groups=df[level])
        else:
            raise ValueError(f"unknown holdout level {level!r}; use 'random' or one of {RANKS}")
        for fold, (_, test) in enumerate(cv):
            out.append(pd.DataFrame({"accession": df.index[test], "level": level, "fold": fold}))
    return pd.concat(out, ignore_index=True).sort_values(["level", "accession"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--usable", required=True)
    ap.add_argument("--levels", required=True)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    g = pd.read_csv(a.genomes, sep="\t").set_index("accession")
    usable = pd.read_csv(a.usable, sep="\t").set_index("accession")["usable"]
    df = g.loc[usable[usable].index.intersection(g.index)].sort_index()
    make_splits(df, a.levels.split(","), a.folds, a.seed).to_csv(a.out, sep="\t", index=False)


if __name__ == "__main__":
    main()

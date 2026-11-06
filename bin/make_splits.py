#!/usr/bin/env python3
"""Benchmark cross-validation folds at nested phylogenetic holdout levels, plus the dev/confirm lockbox.

Usage: make_splits.py --genomes genomes.tsv --levels random,genus,family,order,class --folds 5
                      --seed 42 --dev-fraction 0.5 --out splits.tsv

genomes.tsv: from bench_extract.py (accession, domain..species, group).

level "random": genomes assigned to folds at random (species holdout: each genome is its own GTDB
species). Other levels: whole taxa at that rank are held out together (GroupKFold), so no test genome
shares that taxon with any training genome. Genomes with identical 16S (`group`) always share a
fold: they take the taxonomy of the group's first member for splitting.

partition: families are split once into dev (method development, tuning) and confirm (the frozen
final evaluation, DESIGN §B), stratified by phylum. Every level uses the same partition, so confirm
genomes are never scored during development.

Output: long table accession, level, fold, partition.
"""
import argparse

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def partition(df, dev_fraction, seed):
    """dev/confirm per family, stratified by phylum (each phylum's families shuffled, first share dev)."""
    rng = np.random.default_rng(seed)
    dev = set()
    fams = df[["phylum", "family"]].drop_duplicates().sort_values(["phylum", "family"])
    for _, f in fams.groupby("phylum"):
        names = f.family.to_numpy()[rng.permutation(len(f))]
        dev.update(names[: int(round(len(names) * dev_fraction))])
    return np.where(df.family.isin(dev), "dev", "confirm")


def make_splits(df, levels, folds, seed, dev_fraction=0.5):
    df = df.sort_index()
    if "group" in df:   # identical 16S: split by the group's first member's taxonomy
        first = df.groupby("group").head(1).set_index("group")
        df = df.copy()
        df[RANKS] = first.loc[df.group, RANKS].to_numpy()
    part = partition(df, dev_fraction, seed)
    unit = df["group"] if "group" in df else pd.Series(df.index, index=df.index)
    out = []
    for level in levels:
        if level == "random":
            groups = unit
        elif level in RANKS:
            groups = df[level]
        else:
            raise ValueError(f"unknown holdout level {level!r}; use 'random' or one of {RANKS}")
        # random: shuffled; taxon levels: greedy size balancing (largest taxa first), deterministic
        cv = (GroupKFold(n_splits=folds, shuffle=True, random_state=seed) if level == "random"
              else GroupKFold(n_splits=folds)).split(df, groups=groups)
        for fold, (_, test) in enumerate(cv):
            out.append(pd.DataFrame({"accession": df.index[test], "level": level, "fold": fold,
                                     "partition": part[test]}))
    return pd.concat(out, ignore_index=True).sort_values(["level", "accession"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--levels", required=True)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--dev-fraction", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    g = pd.read_csv(a.genomes, sep="\t").set_index("accession")
    make_splits(g, a.levels.split(","), a.folds, a.seed, a.dev_fraction).to_csv(a.out, sep="\t", index=False)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Filter GTDB species representatives and draw a phylum-stratified sample.

Two modes:
  metadata mode:  select_genomes.py --bac META --arc META [filters] --n N --out genomes.tsv
                  n = 0 keeps every genome that passes the filters (full scale).
  subsample mode: select_genomes.py --from-list genomes.tsv --n 24 --min-archaea 2 --out test.tsv

Allocation (metadata mode, n > 0): each phylum gets min(size, --min-per-phylum) genomes, the rest
is split proportionally to phylum size (largest remainder), then genomes are drawn at random
within each phylum. Everything is seeded, so the output is deterministic.
"""
import argparse
import sys

import numpy as np
import pandas as pd

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
COLS = {
    "accession": "accession",
    "gtdb_taxonomy": "gtdb_taxonomy",
    "checkm2_completeness": "checkm2_completeness",
    "checkm2_contamination": "checkm2_contamination",
    "ssu_length": "ssu_length",
    "ncbi_translation_table": "translation_table",
    "gc_percentage": "gc_percentage",
    "genome_size": "genome_size",
    "gtdb_representative": "gtdb_representative",
}


def load_metadata(paths):
    df = pd.concat(
        [pd.read_csv(p, sep="\t", usecols=list(COLS), dtype=str, low_memory=False) for p in paths],
        ignore_index=True,
    ).rename(columns=COLS)
    for c in ["checkm2_completeness", "checkm2_contamination", "ssu_length", "gc_percentage", "genome_size"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")  # GTDB writes "none" for missing
    df["translation_table"] = pd.to_numeric(df["translation_table"], errors="coerce").fillna(11).astype(int)
    tax = df["gtdb_taxonomy"].str.split(";", expand=True)
    for i, r in enumerate(RANKS):
        df[r] = tax[i].str[3:]  # strip "d__", "p__", ...
    df["ncbi_accession"] = df["accession"].str[3:]  # strip "RS_" / "GB_"
    return df


def filter_genomes(df, min_completeness, max_contamination, min_ssu_len):
    keep = (
        (df["gtdb_representative"] == "t")
        & (df["checkm2_completeness"] >= min_completeness)
        & (df["checkm2_contamination"] <= max_contamination)
        & (df["ssu_length"] >= min_ssu_len)
    )
    return df[keep].copy()


def allocate(sizes, n, floor):
    """Per-group sample sizes: floor first, remainder proportional (largest remainder), capped by size."""
    sizes = pd.Series(sizes)
    alloc = np.minimum(sizes, floor)
    if alloc.sum() > n:
        raise ValueError(f"n={n} is smaller than the phylum floors ({alloc.sum()}); lower --min-per-phylum")
    while alloc.sum() < n:
        room = sizes - alloc
        left = n - alloc.sum()
        share = room / room.sum() * left
        add = np.floor(share).astype(int)
        if add.sum() == 0:  # hand out the last few by largest fractional part
            order = (share - add).sort_values(ascending=False).index
            add[order[:left]] = 1
        alloc = np.minimum(alloc + add, sizes)
    return alloc


def stratified_sample(df, n, floor, seed):
    if n <= 0 or n >= len(df):
        return df
    rng = np.random.default_rng(seed)
    alloc = allocate(df.groupby("phylum").size(), n, floor)
    parts = [
        g.loc[rng.choice(g.index, size=int(alloc[p]), replace=False)]
        for p, g in df.sort_values("accession").groupby("phylum")
        if alloc[p] > 0
    ]
    return pd.concat(parts)


def subsample(df, n, min_archaea, seed):
    rng = np.random.default_rng(seed)
    arc = df[df["domain"] == "Archaea"]
    k = min(min_archaea, len(arc), n)
    pick_arc = rng.choice(arc.index, size=k, replace=False)
    rest = df.index.difference(pick_arc)
    pick_rest = rng.choice(rest, size=min(n - k, len(rest)), replace=False)
    return df.loc[np.concatenate([pick_arc, pick_rest])]


OUT_COLS = ["accession", "ncbi_accession", *RANKS, "translation_table", "checkm2_completeness",
            "checkm2_contamination", "ssu_length", "gc_percentage", "genome_size", "gtdb_taxonomy"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bac", help="bac120_metadata_rXXX.tsv.gz")
    ap.add_argument("--arc", help="ar53_metadata_rXXX.tsv.gz")
    ap.add_argument("--from-list", help="existing genome list to subsample")
    ap.add_argument("--min-completeness", type=float, default=90)
    ap.add_argument("--max-contamination", type=float, default=5)
    ap.add_argument("--min-ssu-len", type=int, default=1400)
    ap.add_argument("--n", type=int, default=0, help="0 = all genomes passing filters")
    ap.add_argument("--min-per-phylum", type=int, default=10)
    ap.add_argument("--min-archaea", type=int, default=2, help="subsample mode only")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    if a.from_list:
        df = pd.read_csv(a.from_list, sep="\t", dtype={"translation_table": int})
        out = subsample(df, a.n, a.min_archaea, a.seed)
    else:
        if not (a.bac and a.arc):
            ap.error("metadata mode needs --bac and --arc")
        df = filter_genomes(load_metadata([a.bac, a.arc]), a.min_completeness, a.max_contamination, a.min_ssu_len)
        print(f"{len(df)} genomes pass filters "
              f"({(df.domain == 'Bacteria').sum()} bacteria, {(df.domain == 'Archaea').sum()} archaea)", file=sys.stderr)
        out = stratified_sample(df, a.n, a.min_per_phylum, a.seed)

    out = out.sort_values("accession")[OUT_COLS]
    out.to_csv(a.out, sep="\t", index=False)
    print(f"wrote {len(out)} genomes from {out.phylum.nunique()} phyla to {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()

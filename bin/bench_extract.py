#!/usr/bin/env python3
"""Benchmark data: turn PICRUSt2's own reference (PICRUSt2-SC) into the benchmark's genome set and truth.

Usage: bench_extract.py --psc psc_dir/ --outdir bench_data/ [--min-ko-genomes 10] [--subsample 0 --seed 42]

psc_dir holds the PICRUSt2 default_files folders `bacteria/` and `archaea/` (copied out of the
PICRUSt2 container): {bac,arc}_ref/ (aligned 16S, tree, HMM, model), ko.txt.gz (KO copy numbers),
16S.txt.gz (16S copy numbers) and *_metadata.csv.gz (GTDB r214 taxonomy).

Both PICRUSt2 and our methods learn from exactly these genomes, sequences and labels (DESIGN §B).
--subsample N > 0 keeps a random N genomes (whole 16S groups), for quick test runs.

Outputs (outdir):
  genomes.tsv     accession, domain, ranks d..s (GTDB r214), ssu_copies, group (identical-16S group:
                  genomes that must share a fold), ssu_len
  ssu.fna         one ungapped full-length 16S per genome (the sequence in PICRUSt2's reference)
  ko_copies.npz   scipy sparse CSR int16, genomes (genomes.tsv order) x KOs (kos.txt order)
  kos.txt         KO ids present in >= --min-ko-genomes genomes
"""
import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
DOMAINS = {"bacteria": "bac_ref", "archaea": "arc_ref"}


def read_aligned(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    seqs[name] = "".join(buf)
    return seqs


def ungap(s):
    return s.replace("-", "").replace(".", "").upper().replace("U", "T")


def groups(accessions, seqs, mates):
    """Union-find over identical 16S sequences and PICRUSt2's '16S cluster' mates."""
    parent = {a: a for a in accessions}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        parent[find(a)] = find(b)

    first = {}
    for a in accessions:
        h = hashlib.sha1(seqs[a].encode()).hexdigest()
        if h in first:
            union(a, first[h])
        else:
            first[h] = a
    for a, m in mates.items():
        for b in m:
            if b in parent:
                union(a, b)
    return {a: find(a) for a in accessions}


def read_domain(psc, domain):
    d = Path(psc) / domain
    ref = DOMAINS[domain]
    aln = read_aligned(d / ref / f"{ref}.fna")
    meta = pd.read_csv(d / f"{domain}_metadata.csv.gz", dtype=str).set_index("accession")
    ko = pd.read_csv(d / "ko.txt.gz", sep="\t", index_col=0)
    ko.columns = [c.removeprefix("ko:") for c in ko.columns]
    ssu_copies = pd.read_csv(d / "16S.txt.gz", sep="\t", index_col=0).iloc[:, 0]
    acc = list(aln)
    assert set(acc) == set(ko.index) == set(ssu_copies.index), f"{domain}: reference files disagree"
    tax = meta.loc[acc, "gtdb_taxonomy"].str.split(";", expand=True)
    tax.columns = RANKS
    tax = tax.apply(lambda c: c.str.slice(3))
    mates = meta.loc[acc, "Other genomes within this 16S cluster"].fillna("")
    df = tax.assign(ssu_copies=ssu_copies.loc[acc].values)
    return df, {a: ungap(aln[a]) for a in acc}, ko.loc[acc], {a: [m for m in v.split(",") if m] for a, v in mates.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--psc", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-ko-genomes", type=int, default=10)
    ap.add_argument("--subsample", type=int, default=0)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    parts = [read_domain(a.psc, d) for d in DOMAINS]
    genomes = pd.concat([p[0] for p in parts])
    seqs = {k: v for p in parts for k, v in p[1].items()}
    mates = {k: v for p in parts for k, v in p[3].items()}
    ko = pd.concat([p[2] for p in parts]).fillna(0)
    assert genomes.index.is_unique, "accession in both domains"

    genomes["group"] = pd.Series(groups(list(genomes.index), seqs, mates))
    genomes["ssu_len"] = [len(seqs[x]) for x in genomes.index]
    genomes.index.name = "accession"
    genomes = genomes.sort_index()
    if 0 < a.subsample < len(genomes):
        pick = np.random.default_rng(a.seed).permutation(genomes.group.unique())
        sizes = genomes.group.value_counts().loc[pick].cumsum()
        genomes = genomes[genomes.group.isin(sizes.index[sizes <= a.subsample])]
    ko = ko.loc[genomes.index]

    keep = ((ko.values >= 1).sum(0) >= a.min_ko_genomes)
    ko = ko.loc[:, keep]
    ko = ko[sorted(ko.columns)]

    genomes.to_csv(out / "genomes.tsv", sep="\t")
    with open(out / "ssu.fna", "w") as f:
        for x in genomes.index:
            f.write(f">{x}\n{seqs[x]}\n")
    sparse.save_npz(out / "ko_copies.npz", sparse.csr_matrix(np.clip(ko.values, 0, 32767).astype(np.int16)))
    (out / "kos.txt").write_text("\n".join(ko.columns) + "\n")
    print(f"{len(genomes)} genomes ({(genomes.domain == 'Archaea').sum()} archaea), {ko.shape[1]} KOs, "
          f"{genomes.group.nunique()} 16S groups")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Patristic distances between usable genomes on the GTDB trees (DESIGN §4.2, distance (b)).

Usage: patristic.py --bac-tree bac120.tree --arc-tree ar53.tree --qc genome_qc.tsv --out patristic.parquet

d(i, j) = K[i, i] + K[j, j] - 2 K[i, j], with K the shared root-to-MRCA length (aim2_prep.kinship).
Bacteria and archaea sit on separate GTDB trees: cross-domain distances are NaN.
shortcut: dense n x n float32 matrix, fine for the pilot (~2k); full scale (~55k) needs a
nearest-neighbour-on-tree algorithm instead (12 GB matrix).
"""
import argparse

import dendropy
import numpy as np
import pandas as pd

from aim2_prep import kinship


def domain_distances(tree_path, keep):
    tree = dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True, rooting="force-rooted")
    names = sorted(keep & {t.label for t in tree.taxon_namespace})
    tree.retain_taxa_with_labels(names)
    tree.suppress_unifurcations()
    K = kinship(tree, names)
    d = np.diag(K)
    return names, (d[:, None] + d[None, :] - 2 * K).astype("float32")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ["--bac-tree", "--arc-tree", "--qc", "--out"]:
        ap.add_argument(a, required=True)
    a = ap.parse_args(argv)
    qc = pd.read_csv(a.qc, sep="\t").set_index("accession")
    keep = set(qc.index[qc.usable])
    blocks = [domain_distances(t, keep) for t in (a.bac_tree, a.arc_tree)]
    names = [n for b in blocks for n in b[0]]
    D = np.full((len(names), len(names)), np.nan, dtype="float32")
    o = 0
    for n, d in blocks:
        D[o:o + len(n), o:o + len(n)] = d
        o += len(n)
    pd.DataFrame(D, index=names, columns=names).to_parquet(a.out)


if __name__ == "__main__":
    main()

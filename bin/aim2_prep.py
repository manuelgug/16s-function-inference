#!/usr/bin/env python3
"""Aim 2 inputs for pyseer: pruned GTDB tree, phylogenetic kinship, allele genotypes, phenotype chunks.

Usage: aim2_prep.py --tree bac120_rXXX.tree --genomes genomes.tsv --qc genome_qc.tsv --features features/
                    --targets targets/ [--domain Bacteria] [--ko-min-prev 0.05] [--ko-max-prev 0.95]
                    [--controls K03170] [--chunks 32] --outdir aim2_input

Kinship K[i, j] = branch length shared from the root to MRCA(i, j); K[i, i] = root-to-tip length.
This is the covariance pyseer's phylogeny_distance --lmm computes (that helper is not shipped in
the pyseer container). Genotypes: the binary allele features, missing = 0, as a pyseer .Rtab.
Phenotypes: KO presence (prevalence filter; control KOs always kept if present in >= 5 genomes) and
module completeness (same filter as Aim 1), split into chunk files for parallel pyseer tasks.
"""
import argparse
from pathlib import Path

import dendropy
import numpy as np
import pandas as pd


def kinship(tree, names):
    """Shared root-to-MRCA branch length for every pair of leaves, in the order of `names`."""
    idx = {n: i for i, n in enumerate(names)}
    K = np.zeros((len(names), len(names)))
    depth = {tree.seed_node: 0.0}
    leaves_under = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            leaves_under[node] = [idx[node.taxon.label]]
        else:
            leaves_under[node] = [i for c in node.child_node_iter() for i in leaves_under[c]]
    for node in tree.preorder_node_iter():  # parents before children, so deeper nodes overwrite
        if node is not tree.seed_node:
            depth[node] = depth[node.parent_node] + (node.edge.length or 0.0)
        L = leaves_under[node]
        K[np.ix_(L, L)] = depth[node]
    return K


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ["--tree", "--genomes", "--qc", "--features", "--targets", "--outdir"]:
        ap.add_argument(a, required=True)
    ap.add_argument("--domain", default="Bacteria")
    ap.add_argument("--ko-min-prev", type=float, default=0.05)
    ap.add_argument("--ko-max-prev", type=float, default=0.95)
    ap.add_argument("--controls", default="K03170")
    ap.add_argument("--chunks", type=int, default=32)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    g = pd.read_csv(a.genomes, sep="\t").set_index("accession")
    qc = pd.read_csv(a.qc, sep="\t").set_index("accession")
    X = pd.read_parquet(Path(a.features) / "features_full.parquet")
    ko = pd.read_parquet(Path(a.targets) / "ko_presence.parquet")
    mod = pd.read_parquet(Path(a.targets) / "module_completeness.parquet")

    tree = dendropy.Tree.get(path=a.tree, schema="newick", preserve_underscores=True, rooting="force-rooted")
    tree_taxa = {t.label for t in tree.taxon_namespace}
    samples = sorted(set(g.index[g.domain == a.domain]) & set(qc.index[qc.usable]) & set(X.index)
                     & set(ko.index) & tree_taxa)
    tree.retain_taxa_with_labels(samples)
    tree.suppress_unifurcations()
    tree.write(path=out / "tree.nwk", schema="newick", suppress_rooting=True)
    pd.DataFrame(kinship(tree, samples), index=samples, columns=samples).to_csv(out / "kinship.tsv", sep="\t")

    geno = X.loc[samples].fillna(0).astype("uint8")
    geno = geno.loc[:, geno.nunique() > 1]
    geno.T.rename_axis("Gene").to_csv(out / "genotypes.Rtab", sep="\t")

    ko, mod = ko.loc[samples], mod.loc[samples]
    prev = ko.mean()
    controls = [c for c in a.controls.split(",") if c in ko and ko[c].sum() >= 5]
    ko = ko.loc[:, prev.between(a.ko_min_prev, a.ko_max_prev) | ko.columns.isin(controls)]
    mod = mod.loc[:, ((mod > 0).mean() >= 0.05) & (mod.var() > 1e-3)]
    for kind, M in [("binary", ko), ("continuous", mod)]:
        for i, cols in enumerate(np.array_split(M.columns, min(a.chunks, M.shape[1]))):
            M[cols].rename_axis("samples").to_csv(out / f"pheno_{kind}_{i:03d}.tsv", sep="\t")
    pd.Series(controls, name="control").to_csv(out / "controls.tsv", sep="\t", index=False)
    print(f"{len(samples)} genomes, {geno.shape[1]} variants, {ko.shape[1]} KOs, {mod.shape[1]} modules, "
          f"controls: {controls}")


if __name__ == "__main__":
    main()

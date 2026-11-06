#!/usr/bin/env python3
"""Build PICRUSt2 custom reference directories for one benchmark fold (training genomes only).

Usage: bench_picrust2_ref.py --psc psc_dir/ --train train_ids.txt --outdir ref/

For each domain, the PICRUSt2-SC reference (aligned 16S, tree, RAxML model, KO and 16S copy tables)
is restricted to the training genomes: alignment rows and table rows are subset, the tree is pruned
to the training tips; constant KO columns are moved to {bac,arc}_constant.tsv. The HMM must be rebuilt from the subset alignment (hmmalign checks that the HMM
came from that MSA); the PICRUSt2 step runs `hmmbuild` (HMMER is in the PICRUSt2 image, not here).
The RAxML model (substitution rates) is kept as shipped: estimated with all genomes, a small bias in
PICRUSt2's favour (DESIGN §B). Held-out genomes are absent from every file; this is asserted.

Outputs: ref/bac_ref/{bac_ref.fna,.tre,.model}, ref/arc_ref/..., ref/{bac,arc}_ko.tsv.gz,
ref/{bac,arc}_constant.tsv, ref/{bac,arc}_16S.tsv.gz
"""
import argparse
import gzip
import shutil
from pathlib import Path

import dendropy
import numpy as np

from bench_extract import DOMAINS


def subset_fasta(src, dst, keep):
    n, write = 0, False
    with open(src) as f, open(dst, "w") as o:
        for line in f:
            if line.startswith(">"):
                write = line[1:].split()[0] in keep
                n += write
            if write:
                o.write(line)
    return n


def subset_table(src, dst, keep, constants=None):
    """Rows of the training genomes, streamed. With `constants` (a path): traits with one value in every training
    genome are removed and listed there (KO, value). castor's max-parsimony HSP fails on any constant trait
    ("invalid 'ncol' value"); every HSP method predicts that same value, which bench_picrust2_pred.py fills in."""
    with gzip.open(src, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in f if line.split("\t", 1)[0] in keep]
    cols = np.arange(len(header))
    if constants is not None:
        vals = np.array([r[1:] for r in rows], dtype=np.float32).reshape(len(rows), len(header) - 1)
        const = (vals.min(0) == vals.max(0)) if rows else np.zeros(len(header) - 1, bool)
        cols = np.concatenate([[0], 1 + np.flatnonzero(~const)])
        with open(constants, "w") as o:
            o.write("ko\tvalue\n" + "".join(f"{header[1 + j]}\t{vals[0, j] if rows else 0:g}\n"
                                            for j in np.flatnonzero(const)))
    with gzip.open(dst, "wt", compresslevel=1) as o:
        for r in [header] + rows:
            o.write("\t".join(r[i] for i in cols) + "\n")
    return len(rows)


def prune_tree(src, dst, keep):
    tree = dendropy.Tree.get(path=src, schema="newick", preserve_underscores=True)
    tree.retain_taxa_with_labels(keep)
    tree.suppress_unifurcations()
    tips = {t.taxon.label for t in tree.leaf_node_iter()}
    tree.write(path=dst, schema="newick", suppress_rooting=True, unquoted_underscores=True)
    return tips


def build(psc, keep, outdir):
    psc, outdir = Path(psc), Path(outdir)
    for domain, ref in DOMAINS.items():
        src, dst = psc / domain / ref, outdir / ref
        dst.mkdir(parents=True, exist_ok=True)
        n_fna = subset_fasta(src / f"{ref}.fna", dst / f"{ref}.fna", keep)
        tips = prune_tree(src / f"{ref}.tre", dst / f"{ref}.tre", keep)
        shutil.copy(src / f"{ref}.model", dst / f"{ref}.model")
        p = ref[:3]
        n_ko = subset_table(psc / domain / "ko.txt.gz", outdir / f"{p}_ko.tsv.gz", keep, outdir / f"{p}_constant.tsv")
        n_16 = subset_table(psc / domain / "16S.txt.gz", outdir / f"{p}_16S.tsv.gz", keep)
        assert tips <= keep and len(tips) == n_fna == n_ko == n_16, f"{domain}: reference files disagree"
        print(f"{domain}: {n_fna} training genomes")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--psc", required=True)
    ap.add_argument("--train", required=True, help="file with one training accession per line")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    build(a.psc, set(Path(a.train).read_text().split()), a.outdir)


if __name__ == "__main__":
    main()

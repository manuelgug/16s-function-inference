#!/usr/bin/env python3
"""Build functional target matrices (genome x target) from eggNOG-mapper annotations.

Usage: build_targets.py --annotations a.emapper.annotations [...] --outdir targets

Outputs (rows = genomes, index "accession"):
  ko_presence.parquet          KO presence/absence (uint8), every KO seen
  module_completeness.parquet  KEGG module completeness in [0, 1] (kegg-pathways-completeness)
  cog_profile.parquet          fraction of COG-annotated proteins per COG category letter;
                               a protein with k letters adds 1/k to each
  annotation_stats.tsv         per-genome counts
Target filtering (prevalence, variance) is done at modelling time on the training set, not here.
Protein ids are "<accession>@<n>" (see process_genomes.py).
"""
import argparse
import io
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


def read_annotations(paths):
    """emapper files start and end with '##' comment lines; the '#query' line is the header."""
    frames = []
    for p in paths:
        with open(p) as f:
            text = "".join(line for line in f if not line.startswith("##"))
        frames.append(pd.read_csv(io.StringIO(text), sep="\t", dtype=str,
                                  usecols=["#query", "COG_category", "KEGG_ko"]))
    df = pd.concat(frames, ignore_index=True)
    df["accession"] = df["#query"].str.rsplit("@", n=1).str[0]
    return df


def ko_sets(df):
    kos = defaultdict(set)
    for acc, ko in zip(df["accession"], df["KEGG_ko"]):
        if isinstance(ko, str) and ko != "-":
            kos[acc].update(k.removeprefix("ko:") for k in ko.split(","))
    return kos


def ko_matrix(kos, genomes):
    rows = [{k: 1 for k in kos.get(g, ())} for g in genomes]
    m = pd.DataFrame(rows, index=pd.Index(genomes, name="accession")).fillna(0).astype("uint8")
    return m.reindex(sorted(m.columns), axis=1)


def cog_matrix(df, genomes):
    acc = defaultdict(lambda: defaultdict(float))
    for g, cats in zip(df["accession"], df["COG_category"]):
        if isinstance(cats, str) and cats != "-":
            for c in cats:
                acc[g][c] += 1 / len(cats)
    m = pd.DataFrame([acc.get(g, {}) for g in genomes], index=pd.Index(genomes, name="accession")).fillna(0.0)
    m = m.reindex(sorted(m.columns), axis=1)
    return m.div(m.sum(axis=1).replace(0, 1), axis=0)


class ModuleScorer:
    """KEGG module completeness exactly as kegg-pathways-completeness computes it, vectorised.

    The library represents a module as a graph whose source-to-sink paths are the alternative ways to
    complete it; each edge is a KO with a weight. Completeness = max over paths of
    (weight of present KOs on the path / total path weight), rounded to 2 decimals in percent; 0 when
    no KO of the module is present. The library re-enumerates every path for every genome (~5 s per
    genome). Here the paths of each module are enumerated once, from the library's own graphs, and
    all genomes are scored with one matrix product per module. Tested equal to `give_completeness`.
    """

    def __init__(self):
        from importlib.resources import files
        import networkx as nx
        from kegg_pathways_completeness.bin.utils import parse_graphs_input
        graphs = parse_graphs_input(files("kegg_pathways_completeness.pathways_data").joinpath("graphs.pkl"))
        self.paths = {}
        for m, (G, _edges, _unnecessary) in graphs.items():
            acc = {}
            for node in nx.topological_sort(G):
                preds = G.pred[node]
                if not preds:
                    acc[node] = [{}]
                    continue
                acc[node] = [{**path, len(path): (e["label"], e["weight"])}
                             for pred in preds for e in preds[pred].values() for path in acc[pred]]
            labels = sorted({lab for path in acc[1] for lab, _ in path.values()})
            col = {lab: i for i, lab in enumerate(labels)}
            W = np.zeros((len(acc[1]), len(labels)))
            for i, path in enumerate(acc[1]):
                for lab, w in path.values():
                    W[i, col[lab]] += w
            self.paths[m] = (labels, W, W.sum(1))

    def matrix(self, presence):
        """presence: DataFrame genomes x KOs (truthy = present) -> DataFrame genomes x modules in [0, 1]."""
        cols = {k: i for i, k in enumerate(presence.columns)}
        P = presence.to_numpy() > 0
        out = {}
        for m, (labels, W, total) in self.paths.items():
            have = np.array([lab in cols for lab in labels])
            if not have.any():
                continue
            sub = P[:, [cols[lab] for lab in np.array(labels)[have]]].astype(float)
            best = (sub @ W[:, have].T / total).max(1)
            pct = np.array([round(b * 100, 2) for b in best])   # the library rounds percent to 2 decimals
            out[m] = np.where(sub.any(1), pct, 0.0) / 100.0
        res = pd.DataFrame(out, index=presence.index).astype("float32")
        res = res.loc[:, res.any()]
        return res.reindex(sorted(res.columns), axis=1)


def module_matrix(kos, genomes):
    return ModuleScorer().matrix(ko_matrix(kos, genomes))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotations", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    df = read_annotations(a.annotations)
    genomes = sorted(df["accession"].unique())
    kos = ko_sets(df)

    ko_matrix(kos, genomes).to_parquet(out / "ko_presence.parquet")
    cog_matrix(df, genomes).to_parquet(out / "cog_profile.parquet")
    module_matrix(kos, genomes).to_parquet(out / "module_completeness.parquet")

    stats = df.groupby("accession").agg(
        n_annotated=("#query", "size"),
        n_with_ko=("KEGG_ko", lambda s: (s != "-").sum()),
        n_with_cog=("COG_category", lambda s: (s != "-").sum()),
    )
    stats["n_distinct_ko"] = pd.Series({g: len(kos.get(g, ())) for g in genomes})
    stats.to_csv(out / "annotation_stats.tsv", sep="\t")


if __name__ == "__main__":
    main()

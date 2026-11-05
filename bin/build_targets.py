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
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

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


def module_matrix(kos, genomes, workdir):
    inp = Path(workdir) / "genome_kos.tsv"
    with open(inp, "w") as f:
        for g in genomes:
            if kos.get(g):
                f.write(g + "\t" + "\t".join(sorted(kos[g])) + "\n")
    out = Path(workdir) / "kpc"
    subprocess.run(["give_completeness", "-i", str(inp), "-m", "-o", str(out)], check=True,
                   stdout=subprocess.DEVNULL)
    long = pd.read_csv(out / "summary.kegg_contigs.tsv", sep="\t", usecols=["contig", "module_accession", "completeness"])
    m = long.pivot(index="contig", columns="module_accession", values="completeness") / 100.0
    m = m.reindex(genomes).fillna(0.0)  # modules not listed have completeness 0
    m.index.name = "accession"
    return m.reindex(sorted(m.columns), axis=1)


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
    with tempfile.TemporaryDirectory() as tmp:
        module_matrix(kos, genomes, tmp).to_parquet(out / "module_completeness.parquet")

    stats = df.groupby("accession").agg(
        n_annotated=("#query", "size"),
        n_with_ko=("KEGG_ko", lambda s: (s != "-").sum()),
        n_with_cog=("COG_category", lambda s: (s != "-").sum()),
    )
    stats["n_distinct_ko"] = pd.Series({g: len(kos.get(g, ())) for g in genomes})
    stats.to_csv(out / "annotation_stats.tsv", sep="\t")


if __name__ == "__main__":
    main()

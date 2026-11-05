#!/usr/bin/env python3
"""Aim 2 summary: calibration, significant 16S allele-function associations, homoplasy, controls.

Usage: aim2_report.py --results *.pyseer.tsv.gz --input aim2_input/ --features features/ --outdir aim2_report

Tests that pyseer flags in `notes` (bad-chisq, high-bse, ...) are excluded from hits.
Significance (DESIGN §5): Bonferroni per phenotype over the number of unique variant patterns
(alpha 0.05), plus Benjamini-Hochberg FDR over all phenotype x variant tests.
Calibration: genomic inflation lambda = median chi2(1) / 0.4549 per phenotype.
Homoplasy: Fitch parsimony changes of each significant allele on the pruned GTDB tree; an
association is interpretable beyond phylogeny only if the allele arose more than once.
"""
import argparse
from pathlib import Path

import dendropy
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import chi2, false_discovery_control  # noqa: E402

COLS = ["phenotype", "variant", "af", "p", "beta", "notes"]
INK2, GRID, SURFACE, C1, C2 = "#52514e", "#e4e3df", "#fcfcfb", "#2a78d6", "#eb6834"


def fitch_changes(tree, geno):
    """Minimum number of state changes per column of a 0/1 genotype table (rows = tree leaves)."""
    idx = {n: i for i, n in enumerate(geno.index)}
    G = geno.to_numpy().astype(np.int8)
    sets, cost = {}, np.zeros(G.shape[1], dtype=int)
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            sets[node] = np.where(G[idx[node.taxon.label]] == 1, 2, 1).astype(np.int8)
            continue
        kids = [sets.pop(c) for c in node.child_node_iter()]
        s = kids[0]
        for k in kids[1:]:  # multifurcations are resolved pairwise
            inter = s & k
            cost += inter == 0
            s = np.where(inter != 0, inter, s | k)
        sets[node] = s
    return pd.Series(cost, index=geno.columns)


def lam(p):
    p = p[(p > 0) & np.isfinite(p)]
    return np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1) if len(p) else np.nan


def md_table(df):
    fmt = lambda v: f"{v:.3g}" if isinstance(v, (float, np.floating)) else str(v)  # noqa: E731
    return "\n".join(["| " + " | ".join(df.columns) + " |", "|" + "---|" * df.shape[1]]
                     + ["| " + " | ".join(fmt(v) for v in r) + " |" for r in df.itertuples(index=False)])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", nargs="+", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--features", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    inp = Path(a.input)
    plt.rcParams.update({"figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
                         "axes.grid": True, "grid.color": GRID, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.edgecolor": GRID, "font.size": 10})

    r = pd.concat(pd.read_csv(f, sep="\t", header=None, names=COLS) for f in a.results if Path(f).stat().st_size)
    r = r.dropna(subset=["p"])
    r["kind"] = np.where(r.phenotype.str.startswith("K"), "KO", "module")
    geno = pd.read_csv(inp / "genotypes.Rtab", sep="\t", index_col=0).T
    n_patterns = geno.T.drop_duplicates().shape[0]
    bonf = 0.05 / n_patterns
    r["q_bh"] = false_discovery_control(r["p"].clip(lower=1e-300), method="bh")
    r["flagged"] = r["notes"].notna() & (r["notes"].astype(str).str.strip() != "")  # pyseer: bad-chisq, high-bse, ...
    r["bonferroni"] = (r["p"] < bonf) & ~r["flagged"]

    per = r.groupby(["phenotype", "kind"]).agg(n_tests=("p", "size"), min_p=("p", "min"),
                                               n_bonferroni=("bonferroni", "sum"),
                                               n_fdr05=("q_bh", lambda q: (q < 0.05).sum()))
    per["lambda"] = r.groupby(["phenotype", "kind"])["p"].apply(lambda p: lam(p.to_numpy()))
    per = per.reset_index()
    per.to_csv(out / "per_phenotype.tsv", sep="\t", index=False)

    hits = r[r.bonferroni].copy()
    meta = pd.read_csv(Path(a.features) / "feature_columns.tsv", sep="\t").set_index("feature")
    hits = hits.join(meta[["ecoli_pos", "state"]], on="variant")
    if len(hits):
        tree = dendropy.Tree.get(path=inp / "tree.nwk", schema="newick", preserve_underscores=True,
                                 rooting="force-rooted")
        changes = fitch_changes(tree, geno[hits.variant.unique()])
        hits["fitch_changes"] = hits.variant.map(changes)
    else:
        hits["fitch_changes"] = pd.Series(dtype=int)
    hits.sort_values("p").to_csv(out / "hits_bonferroni.tsv", sep="\t", index=False)

    # controls: positively associated hit alleles should be G/C-enriched (thermophily)
    controls = pd.read_csv(inp / "controls.tsv", sep="\t")["control"].tolist() if (inp / "controls.tsv").stat().st_size > 8 else []
    ctrl_rows = []
    gc_bg = meta.loc[meta.index.isin(r.variant.unique()), "state"].isin(["G", "C"]).mean()
    for c in controls:
        h = hits[hits.phenotype == c]
        pos = h[h.beta > 0]
        ctrl_rows.append({"control": c, "lambda": per.loc[per.phenotype == c, "lambda"].squeeze(),
                          "n_hits": len(h), "n_positive": len(pos),
                          "gc_fraction_positive_hits": pos.state.isin(["G", "C"]).mean() if len(pos) else np.nan,
                          "gc_fraction_background": gc_bg,
                          "homoplasic_hits(>1 change)": int((h.fitch_changes > 1).sum())})
    ctrl = pd.DataFrame(ctrl_rows)

    # figures
    fig, ax = plt.subplots(figsize=(6, 3.5), constrained_layout=True)
    for kind, color in [("KO", C1), ("module", C2)]:
        ax.hist(per.loc[per.kind == kind, "lambda"].dropna(), bins=40, alpha=0.7, color=color, label=kind)
    ax.axvline(1, color=INK2, lw=1)
    ax.set_xlabel("genomic inflation λ per phenotype (1 = calibrated)")
    ax.set_ylabel("phenotypes")
    ax.legend(frameon=False)
    fig.savefig(out / "fig1_lambda.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 3.5), constrained_layout=True)
    if len(hits):
        c = hits.groupby("ecoli_pos").size()
        ax.bar(c.index, c.values, width=4, color=C1)
    ax.set_xlim(0, 1545)
    ax.set_xlabel("E. coli 16S position")
    ax.set_ylabel("Bonferroni hits (phenotype x allele)")
    fig.savefig(out / "fig2_hits_along_16s.png", dpi=150)
    plt.close(fig)

    homo = hits.drop_duplicates("variant")["fitch_changes"] if len(hits) else pd.Series(dtype=int)
    lines = [
        "# Aim 2 report: 16S alleles associated with function beyond phylogeny", "",
        f"Genomes: {geno.shape[0]} (bacteria). Variants tested: {geno.shape[1]} "
        f"({n_patterns} unique patterns; Bonferroni threshold p < {bonf:.2e}). "
        f"Phenotypes: {(per.kind == 'KO').sum()} KOs, {(per.kind == 'module').sum()} modules.", "",
        "## Calibration", "",
        "![lambda](fig1_lambda.png)", "",
        md_table(per.groupby("kind")["lambda"].describe(percentiles=[.1, .5, .9]).reset_index()
                 [["kind", "count", "10%", "50%", "90%"]]), "",
        "## Significant associations", "",
        f"Bonferroni hits (phenotype × allele): **{len(hits)}** across {hits.phenotype.nunique()} phenotypes "
        f"and {hits.variant.nunique()} alleles. BH-FDR < 0.05 (unflagged): {((r.q_bh < 0.05) & ~r.flagged).sum()}. "
        f"Tests flagged by pyseer (`notes`, e.g. bad-chisq / high-bse) are never hits: "
        f"{r.flagged.mean():.1%} of all tests.", "",
        f"Alleles among hits that arose more than once on the tree (Fitch changes > 1): "
        f"**{int((homo > 1).sum())} of {len(homo)}**. Only these can support association beyond phylogeny.", "",
        "![hits](fig2_hits_along_16s.png)", "",
        "Top 30 hits:", "",
        md_table(hits.sort_values("p").head(30)[["phenotype", "variant", "ecoli_pos", "state", "af", "beta", "p",
                                                 "fitch_changes"]]) if len(hits) else "_none_", "",
        "## Controls", "",
        md_table(ctrl) if len(ctrl) else "_no control phenotype present in >= 5 genomes_", "",
        "Positive control expectation: reverse gyrase (K03170) hits with beta > 0 are G/C-enriched "
        "relative to the background allele composition.", "",
        "Files: `per_phenotype.tsv`, `hits_bonferroni.tsv`.",
    ]
    (out / "aim2_report.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

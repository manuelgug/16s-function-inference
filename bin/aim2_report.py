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
# 16S variable regions, E. coli numbering (Chakravorty et al. 2007, J Microbiol Methods 69:330)
VREGIONS = {"V1": (69, 99), "V2": (137, 242), "V3": (433, 497), "V4": (576, 682), "V5": (822, 879),
            "V6": (986, 1043), "V7": (1117, 1173), "V8": (1243, 1294), "V9": (1435, 1465)}
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
    bins = np.linspace(0, max(2.0, per["lambda"].quantile(0.99)), 41)
    for kind, color in [("KO", C1), ("module", C2)]:
        v = per.loc[per.kind == kind, "lambda"].dropna()
        ax.hist(v, bins=bins, density=True, histtype="step", lw=2, color=color, label=f"{kind} (n={len(v)})")
    ax.axvline(1, color=INK2, lw=1)
    ax.set_xlabel("genomic inflation λ per phenotype (1 = calibrated)")
    ax.set_ylabel("density of phenotypes")
    ax.legend(frameon=False)
    fig.savefig(out / "fig1_lambda.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 3.5), constrained_layout=True)
    for name, (lo, hi) in VREGIONS.items():
        ax.axvspan(lo, hi, color=GRID, alpha=0.6, lw=0, zorder=0)
        ax.text((lo + hi) / 2, 1.01, name, transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                color=INK2, fontsize=9)
    if len(hits):
        c = hits.groupby("ecoli_pos").size()
        ax.bar(c.index, c.values, width=1.5, color=C1, zorder=2)
    ax.set_xlim(0, 1545)
    ax.set_xlabel("E. coli 16S position")
    ax.set_ylabel("Bonferroni hits (phenotype x allele)")
    fig.savefig(out / "fig2_hits_along_16s.png", dpi=150)
    plt.close(fig)

    in_v = lambda pos: pos.apply(lambda x: any(lo <= x <= hi for lo, hi in VREGIONS.values()))  # noqa: E731
    tested_pos = meta.loc[meta.index.isin(r.variant.unique()), "ecoli_pos"].dropna()
    v_tested = in_v(tested_pos).mean() if len(tested_pos) else np.nan
    v_hits = in_v(hits["ecoli_pos"].dropna()).mean() if len(hits) else np.nan
    lam_k = per.loc[per.kind == "KO", "lambda"].dropna()
    cap1 = "\n>\n".join("> " + x for x in [
        "**Figure 1. Is the association test calibrated?** For each phenotype (a KO's presence, or a module's "
        "completeness), the genomic inflation factor λ compares the median association statistic over all tested "
        "alleles with its expectation when nothing is associated. The curves are normalised to the same area, so KOs "
        "and modules can be compared despite different counts.",
        "**How to read it:** λ ≈ 1 means the phylogenetic correction (the linear mixed model with GTDB-tree kinship) "
        "removes the inflation that shared ancestry causes; λ well above 1 would mean leftover population structure "
        "and false positives; λ below 1 means the test is conservative for that phenotype (less power, not false positives).",
        f"**Here:** median λ is {lam_k.median():.2f} for KOs; {lam_k.between(0.8, 1.2).mean():.0%} of KOs have λ "
        "between 0.8 and 1.2."])
    cap2 = "\n>\n".join("> " + x for x in [
        "**Figure 2. Where on the 16S the significant alleles sit.** Each bar counts Bonferroni-significant "
        "phenotype × allele pairs at one alignment position, in *E. coli* 16S numbering. Grey bands mark the nine "
        "variable regions V1-V9; the unshaded stretches are conserved regions.",
        "**How to read it:** tall bars are positions whose alleles track many functions at once. That is typical of "
        "alleles that mark a clade, whose members share many genes: a hotspot is not by itself evidence of a "
        "functional link. Peaks in conserved regions are more surprising than peaks in variable regions.",
        f"**Here:** {v_hits:.0%} of hits fall in variable regions, which hold {v_tested:.0%} of tested alleles."])
    homo = hits.drop_duplicates("variant")["fitch_changes"] if len(hits) else pd.Series(dtype=int)
    lines = [
        "# Aim 2 report: 16S alleles associated with function beyond phylogeny", "",
        f"Genomes: {geno.shape[0]} (bacteria). Variants tested: {geno.shape[1]} "
        f"({n_patterns} unique patterns; Bonferroni threshold p < {bonf:.2e}). "
        f"Phenotypes: {(per.kind == 'KO').sum()} KOs, {(per.kind == 'module').sum()} modules.", "",
        "Each test asks whether genomes that carry one 16S allele (a nucleotide at an aligned position) more often "
        "carry a function, after accounting for shared ancestry. Tool: pyseer linear mixed model; kinship from the "
        "GTDB bacterial tree. One test per phenotype × allele.", "",
        "## Calibration", "",
        "![lambda](fig1_lambda.png)", "", cap1, "",
        md_table(per.groupby("kind")["lambda"].describe(percentiles=[.1, .5, .9]).reset_index()
                 [["kind", "count", "10%", "50%", "90%"]].astype({"count": int})), "",
        "## Significant associations", "",
        f"Bonferroni hits (phenotype × allele): **{len(hits)}** across {hits.phenotype.nunique()} phenotypes "
        f"and {hits.variant.nunique()} alleles. BH-FDR < 0.05 (unflagged): {((r.q_bh < 0.05) & ~r.flagged).sum()}. "
        f"Tests flagged by pyseer (`notes`, e.g. bad-chisq / high-bse) are never hits: "
        f"{r.flagged.mean():.1%} of all tests.", "",
        f"Alleles among hits that arose more than once on the tree (Fitch changes > 1): "
        f"**{int((homo > 1).sum())} of {len(homo)}**. Only these can support association beyond phylogeny.", "",
        "![hits](fig2_hits_along_16s.png)", "", cap2, "",
        "Top 30 hits. `variant` = alignment column_state; `ecoli_pos` = E. coli position; `af` = allele frequency; "
        "`beta` = effect (positive: the allele goes with more of the function); `fitch_changes` = minimum number of "
        "times the allele arose or was lost on the tree.", "",
        md_table(hits.sort_values("p").head(30)[["phenotype", "variant", "ecoli_pos", "state", "af", "beta", "p",
                                                 "fitch_changes"]]) if len(hits) else "_none_", "",
        "## Controls", "",
        "Positive control: reverse gyrase (K03170) is found almost only in hyperthermophiles, whose rRNA is G/C-rich "
        "to stabilise its stems at high temperature. If the method works, alleles positively associated with "
        "K03170 should be G or C more often than tested alleles in general.", "",
        md_table(ctrl) if len(ctrl) else "_no control phenotype present in >= 5 genomes_", "",
        "Positive control expectation: reverse gyrase (K03170) hits with beta > 0 are G/C-enriched "
        "relative to the background allele composition.", "",
        "Files: `per_phenotype.tsv`, `hits_bonferroni.tsv`.",
    ]
    (out / "aim2_report.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Per-genome yield through the pipeline, plus a markdown summary.

Usage: qc_report.py --genomes genomes.tsv --fetch fetch_status.tsv --stats genome_stats.tsv
                    --targets targets/ --features features/ --sina aligned.csv --outdir qc
Writes qc/genome_qc.tsv and qc/qc_report.md.
"""
import argparse
from pathlib import Path

import pandas as pd


def main(argv=None):
    ap = argparse.ArgumentParser()
    for a in ["--genomes", "--fetch", "--stats", "--targets", "--features", "--sina", "--outdir"]:
        ap.add_argument(a, required=True)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    g = pd.read_csv(a.genomes, sep="\t", usecols=["accession", "domain", "phylum"]).set_index("accession")
    fetch = pd.read_csv(a.fetch, sep="\t").set_index("accession")
    stats = pd.read_csv(a.stats, sep="\t").set_index("accession")
    ann = pd.read_csv(Path(a.targets) / "annotation_stats.tsv", sep="\t").set_index("accession")
    feats = pd.read_parquet(Path(a.features) / "features_full.parquet", columns=[])
    sina = pd.read_csv(a.sina).rename(columns={"name": "accession"}).set_index("accession")

    q = g.copy()
    q["fetched"] = (fetch["status"] == "ok").reindex(q.index, fill_value=False)
    q = q.join(stats[["n_proteins", "n_16s", "n_16s_full"]])
    q["has_16s"] = q["n_16s_full"].fillna(0) > 0
    q["annotated"] = ann["n_with_ko"].reindex(q.index).fillna(0) > 0
    q["sina_identity"] = sina["align_ident_slv"].reindex(q.index) if "align_ident_slv" in sina else pd.NA
    q["in_features"] = q.index.isin(feats.index)
    for pcr in sorted(Path(a.features).glob("pcr_*.tsv")):
        arm = pcr.stem.removeprefix("pcr_")
        q[f"pcr_{arm}"] = pd.read_csv(pcr, sep="\t").set_index("accession")["pcr"].reindex(q.index).fillna(0).astype(bool)
    q["usable"] = q["in_features"] & q["annotated"]
    q.to_csv(out / "genome_qc.tsv", sep="\t")

    steps = ["fetched", "has_16s", "annotated", "in_features", "usable"] + [c for c in q if c.startswith("pcr_")]
    by_domain = q.groupby("domain")[steps].sum().T
    by_domain["all"] = q[steps].sum()
    pct = (by_domain["all"] / len(q) * 100).round(1)
    lines = [
        "# QC report", "",
        f"Genomes requested: **{len(q)}** ({', '.join(f'{d}: {n}' for d, n in q.domain.value_counts().items())})", "",
        "| step | " + " | ".join(by_domain.columns) + " | % of requested |",
        "|---|" + "---|" * len(by_domain.columns) + "---|",
        *[f"| {s} | " + " | ".join(str(int(v)) for v in by_domain.loc[s]) + f" | {pct[s]} |" for s in steps],
        "",
        "`usable` = in the 16S feature matrix and has at least one KO annotation; these genomes enter modelling.",
        "`pcr_<arm>` = in-silico PCR success with that primer pair.",
        "",
        f"Design expectation: usable >= 90 % of requested. Observed: **{pct['usable']} %**.",
    ]
    lost = q[~q["usable"]]
    if len(lost):
        lines += ["", f"## Not usable ({len(lost)})", "", "| reason | n |", "|---|---|"]
        reason = lost.apply(lambda r: "download failed" if not r.fetched else "no full-length 16S" if not r.has_16s
                            else "no KO annotation" if not r.annotated else "not in alignment features", axis=1)
        lines += [f"| {k} | {v} |" for k, v in reason.value_counts().items()]
    (out / "qc_report.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

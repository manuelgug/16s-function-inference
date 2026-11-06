#!/usr/bin/env python3
"""Score PICRUSt2 (M2) KO predictions against eggNOG KO presence, on the Aim 1 KO target set.

Usage: picrust2_eval.py --pred KO_predicted.tsv.gz --nsti marker_predicted_and_nsti.tsv.gz --refs reference_genomes.txt
                        --targets targets/ [--ko-min-prev 0.05] [--ko-max-prev 0.95] --out picrust2_metrics.tsv

A KO is predicted present when PICRUSt2's predicted copy number is >= 1. KOs that PICRUSt2 cannot
predict count as absent. Output rows use the Aim 1 genome_metrics schema with model M2_picrust2,
arm full, level external; NSTI is appended as metric "nsti", and "self_in_reference" = 1 when
the genome (GCA/GCF twins count) is one of PICRUSt2's reference genomes (leak).
"""
import argparse
from pathlib import Path

import pandas as pd

from aim1_predict import genome_metrics


def in_reference(genomes, refs_path):
    """True where the genome itself (GCA/GCF twins count) is in PICRUSt2's reference: the leak DESIGN §4.1 warns about."""
    core = lambda s: s.str.replace(r"^(RS_|GB_)?GC[AF]_", "", regex=True).str.split(".").str[0]  # noqa: E731
    refs = set(core(pd.Series(open(refs_path).read().split())))
    return core(pd.Series(genomes)).isin(refs).to_numpy()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ["--pred", "--nsti", "--refs", "--targets", "--out"]:
        ap.add_argument(a, required=True)
    ap.add_argument("--ko-min-prev", type=float, default=0.05)
    ap.add_argument("--ko-max-prev", type=float, default=0.95)
    a = ap.parse_args(argv)

    true = pd.read_parquet(Path(a.targets) / "ko_presence.parquet")
    true = true.loc[:, true.mean().between(a.ko_min_prev, a.ko_max_prev)]
    pred = pd.read_csv(a.pred, sep="\t", index_col=0)
    pred.columns = pred.columns.str.removeprefix("ko:")
    genomes = true.index.intersection(pred.index)
    P = (pred.reindex(index=genomes, columns=true.columns).fillna(0) >= 1).astype(float)
    m = genome_metrics("ko", true.loc[genomes].to_numpy(), P.to_numpy())
    rows = [pd.DataFrame({"accession": genomes, "metric": k, "value": v}) for k, v in m.items()]
    meta = pd.read_csv(a.nsti, sep="\t", index_col=0).reindex(genomes)
    in_ref = in_reference(genomes, a.refs)
    rows.append(pd.DataFrame({"accession": genomes, "metric": "nsti", "value": meta["metadata_NSTI"].to_numpy()}))
    rows.append(pd.DataFrame({"accession": genomes, "metric": "self_in_reference",
                              "value": in_ref.astype(float)}))
    out = pd.concat(rows).assign(fold=-1, model="M2_picrust2", target_type="ko", arm="full", level="external")
    out[["accession", "fold", "model", "target_type", "metric", "value", "arm", "level"]].to_csv(a.out, sep="\t", index=False)


if __name__ == "__main__":
    main()

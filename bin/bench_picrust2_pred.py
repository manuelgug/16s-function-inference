#!/usr/bin/env python3
"""Convert PICRUSt2 predictions of one benchmark task into the shared prediction format.

Usage: bench_picrust2_pred.py --data bench_data/ --task task/ --p2 p2out/ --methods mp,emp_prob,pic,subtree_average
                              --outdir pred/

p2out/ holds, per HSP method m, the KO predictions of each domain's reduced tree
(<m>_ref1.tsv.gz bacteria, <m>_ref2.tsv.gz archaea; every query is in exactly one, its best domain), the
constant traits of each domain's fold reference ({bac,arc}_constant.tsv)
and marker_nsti.tsv.gz (PICRUSt2's combined 16S prediction, NSTI, best domain, closest genome).

PICRUSt2 predicts KO copy numbers. Score = predicted copies; call (as released) = copies >= 1 after
rounding to the nearest integer for the continuous methods (pic, subtree_average), as PICRUSt2's
metagenome step does. KOs PICRUSt2 cannot predict (absent from the fold's reference) and queries it
drops (NSTI > max_nsti) are predicted absent.

Outputs: pred_picrust2_<m>.npz (with copies), nsti.tsv (query, nsti, predicted 16S copies, best_domain,
closest_reference).
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from benchlib import save_pred


def read_domain_tables(p2, method, queries, kos):
    """Each domain's predictions plus its constant traits (dropped from HSP: every method predicts the constant)."""
    parts = []
    for ref, dom in (("ref1", "bac"), ("ref2", "arc")):
        f = Path(p2) / f"{method}_{ref}.tsv.gz"
        if f.exists() and f.stat().st_size > 0:
            t = pd.read_csv(f, sep="\t", index_col=0)
            const = Path(p2) / f"{dom}_constant.tsv"
            if const.exists():
                c = pd.read_csv(const, sep="\t")
                t = t.assign(**dict(zip(c.ko, c.value)))
            t.columns = [c.removeprefix("ko:") for c in t.columns]
            parts.append(t)
    t = pd.concat(parts) if parts else pd.DataFrame()
    assert t.index.is_unique, f"{method}: query predicted in both domains"
    return t.reindex(index=queries, columns=kos).fillna(0).to_numpy(np.float32)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for o in ("--data", "--task", "--p2", "--methods", "--outdir"):
        ap.add_argument(o, required=True)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    kos = (Path(a.data) / "kos.txt").read_text().split()
    queries = pd.read_csv(Path(a.task) / "task.tsv", sep="\t").accession.tolist()

    cols = ["metadata_NSTI", "16S_rRNA_Count", "best_domain", "closest_reference_genome"]
    f = Path(a.p2) / "marker_nsti.tsv.gz"
    nsti = pd.read_csv(f, sep="\t", index_col=0) if f.exists() else pd.DataFrame(columns=cols)
    nsti = nsti.reindex(queries)[cols]
    nsti.columns = ["nsti", "ssu_copies", "best_domain", "closest_reference"]
    nsti.index.name = "query"
    nsti.to_csv(out / "nsti.tsv", sep="\t")

    for m in a.methods.split(","):
        C = read_domain_tables(a.p2, m, queries, kos)
        save_pred(out / f"pred_picrust2_{m}.npz", queries, C, np.rint(C) >= 1, copies=C)


if __name__ == "__main__":
    main()

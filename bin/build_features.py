#!/usr/bin/env python3
"""Turn the SINA alignment into binary allele-at-position features, full-length and per amplicon.

Usage: build_features.py --aln aligned.fasta --ssu ssu.fna --primers primers.tsv --ref-id ecoli_rrsB
                         [--min-col-coverage 0.95] [--min-allele-freq 0.01] [--pcr-error-rate 0.1]
                         --outdir features

Column filter: keep alignment columns where >= min-col-coverage of genomes have a base (not a gap).
States per kept column: A, C, G, T (U is written as T), '-' (gap). IUPAC ambiguity = missing (NaN).
Feature "<ecoli_pos>_<col>_<state>" = 1 if the genome has that state there. Kept if its frequency
is within [f, 1 - f], f = min-allele-freq. ecoli_pos is the E. coli rrsB position at or before the
column (insertions relative to E. coli share its number; the alignment column keeps them unique).

Amplicon arms (one per primer pair in primers.tsv: name, forward, reverse; IUPAC allowed):
in-silico PCR with cutadapt linked primers on each genome's representative 16S. The arm's columns
are those strictly between the median aligned start and end of the amplified inserts (primers
excluded); the arm keeps only genomes whose PCR succeeded, and drops features now invariant.

Outputs: features_full.parquet, features_<arm>.parquet (float32, rows = genomes),
         feature_columns.tsv, pcr_<arm>.tsv, alignment_columns.tsv,
         alignment_kept.parquet (genome x kept column, one character per cell)
"""
import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

STATES = "ACGT-"


def read_alignment(path):
    ids, rows = [], []
    for r in SeqIO.parse(path, "fasta"):
        ids.append(r.id)
        rows.append(np.frombuffer(str(r.seq).upper().replace("U", "T").replace(".", "-").encode(), dtype="S1"))
    if len({len(x) for x in rows}) != 1:
        raise ValueError("alignment rows differ in length")
    return ids, np.vstack(rows)


def ecoli_positions(ref_row):
    """Per column: E. coli position at or before it (0 before the first base)."""
    return np.cumsum(ref_row != b"-")


def ungapped_to_column(row):
    """Index i = column of the i-th base of this aligned row."""
    return np.flatnonzero(row != b"-")


def allele_features(aln, ids, cols, ecoli_pos, min_freq):
    """Binary allele features over the given columns. Returns DataFrame + per-feature metadata."""
    sub = aln[:, cols]
    known = np.isin(sub, [s.encode() for s in STATES])
    data, meta = {}, []
    for j, c in enumerate(cols):
        k = known[:, j]
        n = k.sum()
        if n == 0:
            continue
        for s in STATES:
            hit = (sub[:, j] == s.encode())
            f = hit[k].sum() / n
            if min_freq <= f <= 1 - min_freq and 0 < f < 1:
                name = f"{ecoli_pos[c]}_{c}_{s}"
                data[name] = np.where(k, hit, np.nan).astype("float32")
                meta.append({"feature": name, "column": int(c), "ecoli_pos": int(ecoli_pos[c]),
                             "state": s, "freq": round(float(f), 5)})
    return pd.DataFrame(data, index=pd.Index(ids, name="accession")), pd.DataFrame(meta)


def in_silico_pcr(ssu_fasta, fwd, rev, error_rate, workdir):
    """Return {accession: (insert_start, insert_end)} in ungapped 0-based half-open coordinates."""
    out = Path(workdir) / "insert.fasta"
    rc_rev = str(Seq(rev).reverse_complement())
    subprocess.run(
        ["cutadapt", "-g", f"{fwd};required...{rc_rev};required", "--discard-untrimmed",
         "-e", str(error_rate), "--quiet", "-o", str(out), str(ssu_fasta)],
        check=True,
    )
    seqs = {r.id: str(r.seq).upper() for r in SeqIO.parse(ssu_fasta, "fasta")}
    hits = {}
    for r in SeqIO.parse(out, "fasta"):
        ins = str(r.seq).upper()
        start = seqs[r.id].find(ins)
        if ins and start >= 0:
            hits[r.id] = (start, start + len(ins))
    return hits


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln", required=True)
    ap.add_argument("--ssu", required=True)
    ap.add_argument("--primers", required=True)
    ap.add_argument("--ref-id", default="ecoli_rrsB")
    ap.add_argument("--min-col-coverage", type=float, default=0.95)
    ap.add_argument("--min-allele-freq", type=float, default=0.01)
    ap.add_argument("--pcr-error-rate", type=float, default=0.1)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    ids, aln = read_alignment(a.aln)
    ref = ids.index(a.ref_id)
    ecoli = ecoli_positions(aln[ref])
    keep_rows = [i for i in range(len(ids)) if i != ref]
    ids = [ids[i] for i in keep_rows]
    aln = aln[keep_rows]

    coverage = (aln != b"-").mean(axis=0)
    cols = np.flatnonzero(coverage >= a.min_col_coverage)
    pd.DataFrame({"column": cols, "ecoli_pos": ecoli[cols], "coverage": coverage[cols].round(5)}) \
        .to_csv(out / "alignment_columns.tsv", sep="\t", index=False)
    # kept columns as characters, for 16S identity (k-NN) and per-column analyses (Aim 2)
    pd.DataFrame(aln[:, cols].astype("U1"), index=pd.Index(ids, name="accession"),
                 columns=[str(c) for c in cols]).to_parquet(out / "alignment_kept.parquet")

    full, meta = allele_features(aln, ids, cols, ecoli, a.min_allele_freq)
    full.to_parquet(out / "features_full.parquet")

    primers = pd.read_csv(a.primers, sep="\t", header=None, names=["name", "fwd", "rev"], comment="#")
    row_of = {g: i for i, g in enumerate(ids)}
    with tempfile.TemporaryDirectory() as tmp:
        for p in primers.itertuples():
            hits = in_silico_pcr(a.ssu, p.fwd, p.rev, a.pcr_error_rate, tmp)
            spans = {}
            for g, (s, e) in hits.items():
                if g in row_of:
                    col = ungapped_to_column(aln[row_of[g]])
                    spans[g] = (int(col[s]), int(col[e - 1]))
            pd.DataFrame(
                [{"accession": g, "pcr": int(g in spans), "start_col": spans.get(g, (None,))[0],
                  "end_col": spans.get(g, (None, None))[1]} for g in ids]
            ).to_csv(out / f"pcr_{p.name}.tsv", sep="\t", index=False)
            meta[f"in_{p.name}"] = 0
            if not spans:
                continue
            lo = int(np.median([s for s, _ in spans.values()]))
            hi = int(np.median([e for _, e in spans.values()]))
            in_arm = meta["column"].between(lo, hi, inclusive="neither")
            meta.loc[in_arm, "in_" + p.name] = 1
            arm = full.loc[[g for g in ids if g in spans], meta.loc[in_arm, "feature"]]
            arm = arm.loc[:, arm.nunique(dropna=True) > 1]
            arm.to_parquet(out / f"features_{p.name}.parquet")

    meta.to_csv(out / "feature_columns.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""One benchmark task: training set and query sequences for (level, fold, arm, phase).

Usage: bench_prepare.py --data bench_data/ --splits splits.tsv --level genus --fold 0 --arm V4
                        --primers primers.tsv --pcr-error-rate 0.1 --phase dev --outdir task/

Training genomes = every genome outside the held-out fold (both partitions). Queries = held-out
genomes of the phase's partition (dev or confirm). arm "full": the full-length 16S; otherwise the
in-silico amplicon of the primer pair (genomes whose 16S does not amplify are not queried in that
arm, by any method). Every method, PICRUSt2 included, gets exactly query.fna and the training set.

Leak checks (the task fails if one trips): no query among the training genomes; at a taxon level,
no taxon shared between queries and training genomes.

Outputs: train.txt (training accessions), ref.fna (their full-length 16S: the reference of s16fun and,
in PICRUSt2's own format, of PICRUSt2), train_arm.fna (their arm sequences, for learners trained
per region), query.fna, task.tsv (query accessions with level, fold, arm, phase).
"""
import argparse
import tempfile
from pathlib import Path

import pandas as pd

from benchlib import load_data, write_fasta
from build_features import in_silico_pcr


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for o in ("--data", "--splits", "--level", "--arm", "--primers", "--phase", "--outdir"):
        ap.add_argument(o, required=True)
    ap.add_argument("--fold", type=int, required=True)
    ap.add_argument("--pcr-error-rate", type=float, default=0.1)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    genomes, _, _, seqs = load_data(a.data)
    s = pd.read_csv(a.splits, sep="\t")
    s = s[s.level == a.level]
    test = s[s.fold == a.fold]
    train = sorted(set(genomes.index) - set(test.accession))
    queries = sorted(test.loc[test.partition == a.phase, "accession"])

    assert not set(queries) & set(train), "query genome in the training set"
    if a.level != "random":
        shared = set(genomes.loc[queries, a.level]) & set(genomes.loc[train, a.level])
        assert not shared, f"{len(shared)} {a.level} taxa in both training and query genomes"

    if a.arm == "full":
        arm = seqs
    else:
        pr = pd.read_csv(a.primers, sep="\t", header=None, comment="#").set_index(0)
        fwd, rev = pr.loc[a.arm, 1], pr.loc[a.arm, 2]
        with tempfile.TemporaryDirectory() as tmp:
            write_fasta(Path(tmp) / "all.fna", seqs, sorted(seqs))
            amp = in_silico_pcr(Path(tmp) / "all.fna", fwd, rev, a.pcr_error_rate, tmp)
        arm = {g: seqs[g][s0:e0] for g, (s0, e0) in amp.items()}
    qids = sorted(q for q in queries if q in arm)

    (out / "train.txt").write_text("\n".join(train) + "\n")
    write_fasta(out / "ref.fna", seqs, train)
    write_fasta(out / "query.fna", arm, qids)
    write_fasta(out / "train_arm.fna", arm, [t for t in train if t in arm])
    pd.DataFrame({"accession": qids, "level": a.level, "fold": a.fold, "arm": a.arm, "phase": a.phase}).to_csv(
        out / "task.tsv", sep="\t", index=False)
    print(f"{a.level}/{a.fold}/{a.arm}/{a.phase}: {len(train)} training genomes, {len(qids)} of {len(queries)} queried")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Per-genome processing for one batch: gene calling (pyrodigal) and representative 16S selection.

Usage: process_genomes.py --genomes batch.tsv --fna-dir FNA --gff-dir GFF --prefix B [--min-ssu-len 1400]

Inputs: FNA/<acc>.fna.gz genomes, GFF/<acc>.gff barrnap 0.9 output.
Outputs (all prefixed with --prefix):
  .faa              proteins, ids "<acc>@<n>"
  .ssu.fna          one representative full-length 16S per genome (header = accession)
  .ssu_copies.tsv   every 16S copy found, with length/partial flags and which one was chosen
  .stats.tsv        per-genome counts
The representative copy is the medoid of the genome's full-length, non-partial copies
(minimum summed pairwise distance); ties go to the longest copy.
"""
import argparse
import csv
import gzip
from multiprocessing import Pool
from pathlib import Path

import pyrodigal
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

STATS_COLS = ["accession", "n_contigs", "genome_len", "n_proteins", "gene_mode", "n_16s", "n_16s_full"]
MIN_TRAIN_LEN = 20000  # pyrodigal cannot train single-genome mode on less sequence

_aligner = PairwiseAligner(mode="global", match_score=1, mismatch_score=-1,
                           open_gap_score=-2, extend_gap_score=-1)


def distance(a, b):
    """1 - identities / longer length, from a global alignment."""
    c = _aligner.align(a, b)[0].counts()
    return 1 - c.identities / max(len(a), len(b))


def medoid(seqs):
    """Index of the medoid sequence; ties -> longest, then first."""
    if len(seqs) == 1:
        return 0
    n = len(seqs)
    d = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d[i][j] = d[j][i] = distance(seqs[i], seqs[j])
    return min(range(n), key=lambda i: (round(sum(d[i]), 9), -len(seqs[i]), i))


def parse_barrnap_16s(gff_path):
    hits = []
    for line in open(gff_path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or "Name=16S_rRNA" not in f[8]:
            continue
        hits.append({"contig": f[0], "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                     "partial": "partial" in f[8]})
    return hits


def process(args):
    acc, tt, fna, gff, min_len = args
    with gzip.open(fna, "rt") as fh:
        contigs = {r.id: r.seq for r in SeqIO.parse(fh, "fasta")}

    # genes
    total = sum(len(s) for s in contigs.values())
    if total >= MIN_TRAIN_LEN:
        finder = pyrodigal.GeneFinder(meta=False)
        finder.train(*(bytes(s) for s in contigs.values()), translation_table=tt)
    else:
        finder = pyrodigal.GeneFinder(meta=True)
    proteins, k = [], 0
    for s in contigs.values():
        for g in finder.find_genes(bytes(s)):
            k += 1
            proteins.append(f">{acc}@{k}\n{g.translate(include_stop=False)}\n")

    # 16S copies
    copies = []
    for i, h in enumerate(parse_barrnap_16s(gff) if Path(gff).exists() else []):
        seq = contigs[h["contig"]][h["start"] - 1:h["end"]]
        seq = str(seq.reverse_complement() if h["strand"] == "-" else seq).upper()
        copies.append({**h, "copy": i + 1, "length": len(seq), "seq": seq})
    full = [c for c in copies if c["length"] >= min_len and not c["partial"]]
    chosen = full[medoid([c["seq"] for c in full])]["copy"] if full else None

    stats = {"accession": acc, "n_contigs": len(contigs), "genome_len": total, "n_proteins": k,
             "gene_mode": "single" if total >= MIN_TRAIN_LEN else "meta",
             "n_16s": len(copies), "n_16s_full": len(full)}
    return acc, proteins, copies, chosen, stats


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--fna-dir", required=True)
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--min-ssu-len", type=int, default=1400)
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args(argv)

    jobs = []
    with open(a.genomes) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            fna = Path(a.fna_dir) / f"{r['accession']}.fna.gz"
            if fna.exists():  # missing = failed download, already reported by fetch
                jobs.append((r["accession"], int(r["translation_table"]), fna,
                             Path(a.gff_dir) / f"{r['accession']}.gff", a.min_ssu_len))

    with Pool(a.threads) as pool, \
            open(f"{a.prefix}.faa", "w") as faa, open(f"{a.prefix}.ssu.fna", "w") as ssu, \
            open(f"{a.prefix}.ssu_copies.tsv", "w") as cp, open(f"{a.prefix}.stats.tsv", "w") as st:
        cp.write("accession\tcopy\tcontig\tstart\tend\tstrand\tlength\tpartial\tchosen\n")
        st.write("\t".join(STATS_COLS) + "\n")
        for acc, proteins, copies, chosen, stats in pool.imap(process, jobs):
            faa.writelines(proteins)
            for c in copies:
                cp.write(f"{acc}\t{c['copy']}\t{c['contig']}\t{c['start']}\t{c['end']}\t{c['strand']}\t"
                         f"{c['length']}\t{int(c['partial'])}\t{int(c['copy'] == chosen)}\n")
                if c["copy"] == chosen:
                    ssu.write(f">{acc}\n{c['seq']}\n")
            st.write("\t".join(str(stats[h]) for h in STATS_COLS) + "\n")


if __name__ == "__main__":
    main()

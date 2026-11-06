#!/usr/bin/env python3
"""k-mer neural network baseline (embeRNA family) on one benchmark task.

Usage: bench_mlp.py --data bench_data/ --task task/ [--hidden 512,512] [--max-iter 60] [--seed 42]
                    --outdir pred/

Re-implementation of the embeRNA design (bioRxiv 2026.04.16.718937) on this benchmark's labels:
1- to 5-mer frequencies of the 16S (1,364 features, each k normalised to sum 1) -> two hidden layers
of 512 ReLU units -> one sigmoid per KO, binary cross-entropy, Adam (lr 1e-3, batch 1,024), stopped
when the training loss improves < 1e-4 for 10 epochs or after --max-iter epochs (scikit-learn's early
stopping is not used: it scores multi-label validation by exact-match accuracy, which is always 0 with
thousands of KOs, and stops training after 10 epochs). Trained on the training genomes' sequences of the task's
arm (full length or the amplicon), as embeRNA trains one model per region. KOs absent from every
training genome are predicted absent.

Output: pred_kmer_mlp.npz (score = probability, call = probability >= 0.5).
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.neural_network import MLPClassifier

from benchlib import load_data, read_fasta, save_pred


def kmer_features(seqs, kmax=5):
    """Rows: 1..kmax-mer frequencies (each k normalised separately); k-mers with non-ACGT are skipped.
    Column order per k: lexicographic over ACGT (A = 0 ... T = 3, first base most significant)."""
    lut = np.full(256, -1, np.int64)
    lut[np.frombuffer(b"ACGT", np.uint8)] = np.arange(4)
    blocks = [np.zeros((len(seqs), 4 ** k), np.float32) for k in range(1, kmax + 1)]
    for r, s in enumerate(seqs):
        x = lut[np.frombuffer(s.upper().encode(), np.uint8)]
        for k in range(1, kmax + 1):
            if len(x) < k:
                continue
            w = np.lib.stride_tricks.sliding_window_view(x, k)
            w = w[(w >= 0).all(1)]
            blocks[k - 1][r] = np.bincount(w @ (4 ** np.arange(k - 1, -1, -1)), minlength=4 ** k)
    return np.hstack([b / np.maximum(b.sum(1, keepdims=True), 1) for b in blocks])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True)
    ap.add_argument("--task", required=True)
    ap.add_argument("--hidden", default="512,512")
    ap.add_argument("--max-iter", type=int, default=60)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    out, task = Path(a.outdir), Path(a.task)
    out.mkdir(parents=True, exist_ok=True)

    genomes, kos, copies, _ = load_data(a.data)
    row = {g: i for i, g in enumerate(genomes.index)}
    train = read_fasta(task / "train_arm.fna")
    queries = pd.read_csv(task / "task.tsv", sep="\t").accession.tolist()
    qseqs = read_fasta(task / "query.fna")
    tids = sorted(train)
    Y = (copies[[row[t] for t in tids]] > 0).toarray()
    present = Y.any(0)

    net = MLPClassifier(hidden_layer_sizes=tuple(int(h) for h in a.hidden.split(",")), activation="relu",
                        solver="adam", learning_rate_init=1e-3, batch_size=1024, max_iter=a.max_iter,
                        early_stopping=False, tol=1e-4, n_iter_no_change=10, random_state=a.seed)
    net.fit(kmer_features([train[t] for t in tids]), Y[:, present])
    P = np.zeros((len(queries), len(kos)), np.float32)
    if queries:
        P[:, present] = net.predict_proba(kmer_features([qseqs[q] for q in queries]))
    save_pred(out / "pred_kmer_mlp.npz", queries, P, P >= 0.5)
    print(f"{len(tids)} training sequences, {present.sum()} KOs, {net.n_iter_} epochs")


if __name__ == "__main__":
    main()

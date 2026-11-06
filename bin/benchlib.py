"""Shared helpers for the benchmark scripts (bench_*.py): data loading, prediction files, metrics."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


# ---------- data ----------

def load_data(d):
    """bench_extract.py outputs -> (genomes DataFrame, KO list, copies CSR, ungapped 16S dict)."""
    d = Path(d)
    g = pd.read_csv(d / "genomes.tsv", sep="\t", index_col=0)
    kos = (d / "kos.txt").read_text().split()
    copies = sparse.load_npz(d / "ko_copies.npz").tocsr()
    return g, kos, copies, read_fasta(d / "ssu.fna")


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(path, seqs, ids):
    with open(path, "w") as f:
        for i in ids:
            f.write(f">{i}\n{seqs[i]}\n")


# ---------- prediction files ----------
# One file per (task, method): queries (str), score (float16, queries x KOs, higher = more likely
# present), call (bool, packed), optional copies (float16). KO order = kos.txt.

def save_pred(path, queries, score, call, copies=None):
    arrays = dict(queries=np.asarray(queries, dtype=str), score=np.asarray(score, np.float16),
                  call=np.packbits(np.asarray(call, bool), axis=1), n_ko=np.int64(np.shape(call)[1]))
    if copies is not None:
        arrays["copies"] = np.asarray(copies, np.float16)
    np.savez_compressed(path, **arrays)


def load_pred(path):
    z = np.load(path)
    out = dict(queries=z["queries"].tolist(), score=z["score"].astype(np.float32),
               call=np.unpackbits(z["call"], axis=1, count=int(z["n_ko"])).astype(bool))
    if "copies" in z:
        out["copies"] = z["copies"].astype(np.float32)
    return out


# ---------- metrics ----------

def ko_scores(Y, P):
    """Per-row KO-set scores for boolean truth Y and calls P (rows = genomes)."""
    tp = (Y & P).sum(1)
    fp = (~Y & P).sum(1)
    fn = (Y & ~P).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        prec = np.where(tp + fp > 0, tp / (tp + fp), 0.0)
        rec = np.where(tp + fn > 0, tp / (tp + fn), 0.0)
        f1 = np.where(2 * tp + fp + fn > 0, 2 * tp / (2 * tp + fp + fn), 1.0)
        jac = np.where(tp + fp + fn > 0, tp / (tp + fp + fn), 1.0)
    return dict(ko_f1=f1, ko_precision=prec, ko_recall=rec, ko_jaccard=jac, ko_tp=tp, ko_fp=fp, ko_fn=fn,
                ko_n_predicted=tp + fp, ko_n_true=tp + fn)


def cluster_bootstrap(diff, clusters, rng, n_boot=2000):
    """Mean of per-genome paired differences with a 95 % CI from resampling whole clusters (e.g. genera):
    genomes of one cluster are not independent, so they are resampled together."""
    diff = np.asarray(diff, float)
    codes, uniq = pd.factorize(np.asarray(clusters))
    sums = np.bincount(codes, weights=diff, minlength=len(uniq))
    counts = np.bincount(codes, minlength=len(uniq))
    draw = rng.integers(0, len(uniq), size=(n_boot, len(uniq)))
    boot = sums[draw].sum(1) / counts[draw].sum(1)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return diff.mean(), lo, hi


# ---------- per-KO AUROC under bootstrap weights ----------

class WeightedAUROC:
    """Per-column AUROC of scores S against truth Y (genomes x KOs) for any genome weights, exact with ties
    (a tie counts 1/2). Sorting and tie groups are computed once; each weighting then costs two bincounts
    over the flattened matrix, so a cluster bootstrap over thousands of KOs is affordable."""

    def __init__(self, S, Y):
        S = np.asarray(S, np.float32)
        Y = np.asarray(Y, bool)
        n, k = S.shape
        order = np.argsort(S, axis=0, kind="stable")
        Ss = np.take_along_axis(S, order, axis=0)
        new = np.vstack([np.ones((1, k), bool), Ss[1:] != Ss[:-1]])      # start of a tie group
        gid = np.cumsum(new, axis=0) - 1                                   # group index within column
        ngroups = gid[-1] + 1
        offset = np.concatenate([[0], np.cumsum(ngroups)[:-1]])
        self.flat_gid = (gid + offset[None, :]).ravel(order="F")           # column-major: one column after another
        self.rows = order.ravel(order="F")
        self.ypos = np.take_along_axis(Y, order, axis=0).ravel(order="F")
        self.ngroups, self.offset, self.k, self.n = ngroups, offset, k, n
        self.col_of_group = np.repeat(np.arange(k), ngroups)

    def __call__(self, w):
        w = np.asarray(w, np.float64)[self.rows]
        G = int(self.ngroups.sum())
        wpos = np.bincount(self.flat_gid, weights=w * self.ypos, minlength=G)
        wneg = np.bincount(self.flat_gid, weights=w * ~self.ypos, minlength=G)
        cneg = np.cumsum(wneg)
        start = cneg[self.offset] - wneg[self.offset]                      # negatives before each column
        neg_before = cneg - wneg - start[self.col_of_group]
        num = np.bincount(self.col_of_group, weights=wpos * (neg_before + 0.5 * wneg), minlength=self.k)
        P = np.bincount(self.col_of_group, weights=wpos, minlength=self.k)
        N = np.bincount(self.col_of_group, weights=wneg, minlength=self.k)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where((P > 0) & (N > 0), num / (P * N), np.nan)

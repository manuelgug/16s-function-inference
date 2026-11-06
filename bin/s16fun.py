#!/usr/bin/env python3
"""s16fun: predict KEGG functions (KOs, modules) of prokaryotes from 16S rRNA sequences.

  s16fun.py build-ref --ssu ssu.fna --genomes genomes.tsv --ko-copies ko_copies.npz --kos kos.txt
                      [--model stacker.txt] --out ref/
  s16fun.py predict --ref ref/ --seqs asv.fna [--table counts.tsv] --out out/ [--threads 8]

Method (DESIGN §B):
  1. search: vsearch global alignment of each query against the reference 16S; the top K hits with
     their identity (no multiple alignment, no all-pairs matrix; amplicons work as they are).
  2. evidence per (query, KO): what the nearest genome has (1-NN), identity-weighted averages over
     the 5 and 25 nearest genomes, the KO's prevalence in the nearest genome's genus, family, order,
     class and phylum, the KO's global prevalence, the KO's phylogenetic conservation, and how close
     the nearest genomes are (top identity, mean of the top 5).
  3. a gradient-boosted stacker turns the evidence into a calibrated probability. It learns how much
     to trust close neighbours versus clade averages as the nearest reference gets farther away.
  4. calls: per query, the KO set that maximises expected F1 under the probabilities.
  5. 16S copy number (for community profiles): a second learned model on the copy numbers and identities of
     the 10 nearest genomes and clade means; community KO abundance = sum over ASVs of reads / 16S copies x
     expected KO copies (probability x copies when present).
Without models, step 3 falls back to the 5-NN average and step 5 to the 5-NN weighted median.

Reference dir: ssu.fna, genomes.tsv (accession, domain..species, ssu_copies), kos.txt,
ko_copies.npz (genomes x KOs, copy numbers), [stacker.txt], [ssu_model.txt].

predict outputs (out/): ko_probability.tsv.gz, ko_calls.tsv.gz (0/1), ko_copies.tsv.gz (expected copies),
module_completeness.tsv.gz, novelty.tsv (top identity, nearest genome, predicted 16S copies),
and with --table: community_ko.tsv.gz (16S-copy-normalised, copy-number weighted KO abundance per sample).
"""
import argparse
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

PRIOR_RANKS = ["genus", "family", "order", "class", "phylum"]
FEATURES = ["knn1", "knn5", "knn25"] + [f"prior_{r}" for r in PRIOR_RANKS] + [
    "prevalence", "conservation", "top_id", "mean_id5", "id_gap"]
K = 25


# ---------- reference ----------

class Reference:
    def __init__(self, genomes, ssu_fna, kos, copies):
        self.genomes = genomes                      # DataFrame indexed by accession, row order = matrices
        self.ssu_fna = Path(ssu_fna) if ssu_fna else None
        self.kos = list(kos)
        self.copies = sparse.csr_matrix(copies)
        self.Y = (self.copies > 0).astype(np.float32).tocsr()
        self.row = {a: i for i, a in enumerate(genomes.index)}
        self._clade = {}
        self.prevalence = np.asarray(self.Y.mean(0)).ravel().astype(np.float32)
        self.conservation = self._conservation()

    @classmethod
    def load(cls, d):
        d = Path(d)
        g = pd.read_csv(d / "genomes.tsv", sep="\t", index_col=0)
        return cls(g, d / "ssu.fna", (d / "kos.txt").read_text().split(), sparse.load_npz(d / "ko_copies.npz"))

    def subset(self, accessions, ssu_fna):
        """Reference restricted to `accessions` (benchmark folds); ssu_fna must hold exactly them."""
        idx = [self.row[a] for a in accessions]
        ids = [l[1:].split()[0] for l in open(ssu_fna) if l.startswith(">")]
        assert sorted(ids) == sorted(accessions), "reference FASTA and training genomes disagree"
        return Reference(self.genomes.iloc[idx], ssu_fna, self.kos, self.copies[idx])

    def _conservation(self):
        """Share of a KO's variance that lies between genera (1 = constant within genera)."""
        g = pd.factorize(self.genomes["genus"])[0]
        M = sparse.csr_matrix((np.ones(len(g), np.float32), (g, np.arange(len(g)))))
        size = np.asarray(M.sum(1)).ravel()
        mean_g = (M @ self.Y).toarray() / size[:, None]
        within = (size[:, None] * mean_g * (1 - mean_g)).sum(0) / len(g)
        total = self.prevalence * (1 - self.prevalence)
        return np.where(total > 0, 1 - within / np.maximum(total, 1e-12), 1.0).astype(np.float32)

    def clade_prevalence(self, rank, taxa):
        """KO prevalence in each requested taxon at `rank` (rows follow `taxa`; taxa from this reference)."""
        if rank not in self._clade:
            codes, uniq = pd.factorize(self.genomes[rank])
            M = sparse.csr_matrix((np.ones(len(codes), np.float32), (codes, np.arange(len(codes)))))
            prev = (M @ self.Y).toarray() / np.asarray(M.sum(1))
            self._clade[rank] = (pd.Index(uniq), prev.astype(np.float32))
        uniq, prev = self._clade[rank]
        return prev[uniq.get_indexer(taxa)]


# ---------- search ----------

def search(query_fna, db_fna, threads, workdir, k=K):
    """Top-k vsearch hits per query: DataFrame query, target, identity (fraction), sorted by identity."""
    out = Path(workdir) / "hits.tsv"
    subprocess.run(["vsearch", "--usearch_global", str(query_fna), "--db", str(db_fna), "--id", "0.5",
                    "--maxaccepts", str(k), "--maxrejects", str(max(256, 8 * k)), "--strand", "plus",
                    "--userout", str(out), "--userfields", "query+target+id", "--threads", str(threads),
                    "--quiet"], check=True)
    h = pd.read_csv(out, sep="\t", names=["query", "target", "identity"], dtype={"query": str, "target": str})
    h["identity"] /= 100.0
    return h.sort_values(["query", "identity", "target"], ascending=[True, False, True], kind="stable")


# ---------- evidence ----------

def knn_average(Y, nbr, w, k):
    """Weighted average of neighbour rows; nbr/w: queries x K arrays (nbr -1 = no hit)."""
    n, kk = nbr[:, :k], w[:, :k]
    out = np.zeros((nbr.shape[0], Y.shape[1]), np.float32)
    tot = kk.sum(1)
    for j in range(n.shape[1]):
        ok = n[:, j] >= 0
        if ok.any():
            out[ok] += kk[ok, j, None] * Y[n[ok, j]].toarray()
    return out / np.where(tot > 0, tot, 1)[:, None]


def neighbours(hits, queries, ref, k=K):
    nbr = np.full((len(queries), k), -1)
    ids = np.zeros((len(queries), k), np.float32)
    qpos = {q: i for i, q in enumerate(queries)}
    for q, grp in hits.groupby("query", sort=False):
        if q not in qpos:
            continue
        t = grp.head(k)
        nbr[qpos[q], :len(t)] = [ref.row[x] for x in t.target]
        ids[qpos[q], :len(t)] = t.identity.to_numpy()
    w = np.where(nbr >= 0, 1.0 / (1.0 - ids + 0.001), 0.0).astype(np.float32)
    return nbr, ids, w


def evidence(hits, queries, ref):
    """dict of queries x KOs evidence matrices (FEATURES) plus per-query novelty."""
    nbr, ids, w = neighbours(hits, queries, ref)
    hit = nbr[:, 0] >= 0
    top = np.where(hit, nbr[:, 0], 0)
    ev = {"knn1": np.where(hit[:, None], ref.Y[top].toarray(), ref.prevalence[None, :]).astype(np.float32),
          "knn5": knn_average(ref.Y, nbr, w, 5), "knn25": knn_average(ref.Y, nbr, w, K)}
    for r in PRIOR_RANKS:
        ev[f"prior_{r}"] = np.where(hit[:, None], ref.clade_prevalence(r, ref.genomes[r].to_numpy()[top]),
                                    ref.prevalence[None, :])
    nq, nk = len(queries), len(ref.kos)
    ev["prevalence"] = np.broadcast_to(ref.prevalence, (nq, nk))
    ev["conservation"] = np.broadcast_to(ref.conservation, (nq, nk))
    n5 = np.maximum((nbr[:, :5] >= 0).sum(1), 1)
    novelty = pd.DataFrame({"top_id": ids[:, 0], "mean_id5": ids[:, :5].sum(1) / n5,
                            "id_gap": ids[:, 0] - ids[:, min(4, K - 1)],
                            "nearest": [ref.genomes.index[i] if i >= 0 else "" for i in nbr[:, 0]],
                            "ssu_copies": weighted_median(ref.genomes.ssu_copies.to_numpy(float), nbr, w, 5)},
                           index=pd.Index(queries, name="query"))
    for c in ["top_id", "mean_id5", "id_gap"]:
        ev[c] = np.broadcast_to(novelty[c].to_numpy(np.float32)[:, None], (nq, nk))
    ev["_cond_copies"] = conditional_copies(ref, nbr, w)
    ev["_ssu_X"] = ssu_features(ref, nbr, ids, w)
    return ev, novelty


def weighted_median(values, nbr, w, k):
    """Identity-weighted median of a reference trait over the k nearest genomes (1 when there is no hit)."""
    out = np.ones(len(nbr))
    for i in range(len(nbr)):
        ok = nbr[i, :k] >= 0
        if ok.any():
            v, ww = values[nbr[i, :k][ok]], w[i, :k][ok]
            o = np.argsort(v)
            out[i] = v[o][np.searchsorted(np.cumsum(ww[o]), ww.sum() / 2)]
    return out


SSU_K = 10


def ssu_features(ref, nbr, ids, w):
    """Features of the 16S copy-number model: log copies and identity of the 10 nearest genomes, the log
    geometric mean of copies in the nearest genome's genus and family, and the 5-NN weighted median."""
    c = np.log(ref.genomes.ssu_copies.to_numpy(float))
    k = nbr[:, :SSU_K]
    lc = np.where(k >= 0, c[np.maximum(k, 0)], np.nan)
    top = np.maximum(nbr[:, 0], 0)
    clade = [ref.genomes.assign(lc=c).groupby(r).lc.mean().reindex(ref.genomes[r].to_numpy()[top]).to_numpy()
             for r in ("genus", "family")]
    med = np.log(weighted_median(np.exp(c), nbr, w, 5))
    return np.column_stack([lc, ids[:, :SSU_K]] + clade + [med]).astype(np.float32)


def ssu_copies(ev, novelty, model=None):
    """16S copy numbers: the learned model when given (trained on log copies), else the 5-NN weighted median."""
    if model is None:
        return novelty.ssu_copies.to_numpy()
    return np.maximum(np.exp(model.predict(ev["_ssu_X"])), 1.0)


def conditional_copies(ref, nbr, w):
    """Expected copy number of a KO given that it is present: weighted mean over the 25 nearest
    genomes that carry it (1 when none does)."""
    num = knn_average(ref.copies.astype(np.float32).tocsr(), nbr, w, K)
    den = knn_average(ref.Y, nbr, w, K)
    return np.where(den > 0, num / np.maximum(den, 1e-12), 1.0).astype(np.float32)


def feature_rows(ev, rows=None):
    """Stack evidence into a (queries*KOs) x len(FEATURES) matrix (row-major: query, then KO)."""
    cols = [ev[f] if rows is None else ev[f][rows] for f in FEATURES]
    return np.stack([np.asarray(c, np.float32).ravel() for c in cols], axis=1)


# ---------- probabilities and calls ----------

def probabilities(ev, stacker=None, chunk=2_000_000):
    if stacker is None:
        return ev["knn5"]
    nq, nk = ev["knn5"].shape
    per = max(1, chunk // nk)
    out = np.empty((nq, nk), np.float32)
    for s in range(0, nq, per):
        rows = slice(s, min(nq, s + per))
        out[rows] = stacker.predict(feature_rows(ev, rows)).reshape(-1, nk)
    return out


def expected_f1_calls(P):
    """Per row, call the top-n KOs with n maximising the plug-in expected F1 2*sum(p[:n]) / (n + sum(p))."""
    order = np.argsort(-P, axis=1, kind="stable")
    ps = np.take_along_axis(P, order, axis=1)
    csum = np.cumsum(ps, axis=1)
    n = np.arange(1, P.shape[1] + 1)[None, :]
    ef1 = 2 * csum / (n + P.sum(1, keepdims=True))
    best = ef1.argmax(1)
    calls = np.zeros_like(P, dtype=bool)
    keep = n <= (best[:, None] + 1)
    np.put_along_axis(calls, order, keep, axis=1)
    return calls


def load_stacker(path):
    if path is None or not Path(path).exists():
        return None
    import lightgbm as lgb
    return lgb.Booster(model_file=str(path))


def predict(ref, query_fna, threads, workdir, stacker=None, ssu_model=None):
    """Full prediction for the sequences in query_fna: probabilities, calls, copies, novelty."""
    queries = [l[1:].split()[0] for l in open(query_fna) if l.startswith(">")]
    hits = search(query_fna, ref.ssu_fna, threads, workdir)
    ev, novelty = evidence(hits, queries, ref)
    novelty["ssu_copies"] = ssu_copies(ev, novelty, ssu_model)
    P = probabilities(ev, stacker)
    calls = expected_f1_calls(P)
    return dict(queries=queries, P=P, calls=calls, copies=P * ev["_cond_copies"], novelty=novelty, ev=ev, hits=hits)


# ---------- CLI ----------

def build_ref(a):
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    g = pd.read_csv(a.genomes, sep="\t", index_col=0)
    seqs = {}
    name = None
    for line in open(a.ssu):
        if line.startswith(">"):
            name = line[1:].split()[0]
        else:
            seqs[name] = seqs.get(name, "") + line.strip()
    assert list(g.index) == sorted(g.index) and set(seqs) == set(g.index), "genomes.tsv and ssu.fna disagree"
    with open(out / "ssu.fna", "w") as f:
        for x in g.index:
            f.write(f">{x}\n{seqs[x]}\n")
    g.to_csv(out / "genomes.tsv", sep="\t")
    sparse.save_npz(out / "ko_copies.npz", sparse.load_npz(a.ko_copies))
    Path(out / "kos.txt").write_text(Path(a.kos).read_text())
    if a.model:
        Path(out / "stacker.txt").write_text(Path(a.model).read_text())
    if a.ssu_model:
        Path(out / "ssu_model.txt").write_text(Path(a.ssu_model).read_text())


def run_predict(a):
    from build_targets import ModuleScorer
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ref = Reference.load(a.ref)
    with tempfile.TemporaryDirectory() as tmp:
        r = predict(ref, a.seqs, a.threads, tmp, load_stacker(Path(a.ref) / "stacker.txt"),
                    load_stacker(Path(a.ref) / "ssu_model.txt"))
    idx = pd.Index(r["queries"], name="query")
    frame = lambda m, fmt: pd.DataFrame(m, index=idx, columns=ref.kos)
    frame(r["P"].round(4), None).to_csv(out / "ko_probability.tsv.gz", sep="\t")
    calls = frame(r["calls"].astype(np.uint8), None)
    calls.to_csv(out / "ko_calls.tsv.gz", sep="\t")
    copies = frame(r["copies"], None)
    copies.round(3).to_csv(out / "ko_copies.tsv.gz", sep="\t")
    ModuleScorer().matrix(calls).round(4).to_csv(out / "module_completeness.tsv.gz", sep="\t")
    r["novelty"].to_csv(out / "novelty.tsv", sep="\t")
    if a.table:
        counts = pd.read_csv(a.table, sep="\t", index_col=0).reindex(idx).fillna(0)
        cells = counts.div(r["novelty"]["ssu_copies"].clip(lower=1), axis=0)   # 16S copy normalisation
        (cells.T @ copies).T.round(3).to_csv(out / "community_ko.tsv.gz", sep="\t")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build-ref")
    for o in ("--ssu", "--genomes", "--ko-copies", "--kos", "--out"):
        b.add_argument(o, required=True)
    b.add_argument("--model", help="KO stacker (LightGBM text model)")
    b.add_argument("--ssu-model", help="16S copy-number model (LightGBM text model)")
    p = sub.add_parser("predict")
    p.add_argument("--ref", required=True)
    p.add_argument("--seqs", required=True)
    p.add_argument("--table")
    p.add_argument("--out", required=True)
    p.add_argument("--threads", type=int, default=os.cpu_count())
    a = ap.parse_args(argv)
    {"build-ref": build_ref, "predict": run_predict}[a.cmd](a)


if __name__ == "__main__":
    main()

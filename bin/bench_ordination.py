"""Ordination of predicted against true functional profiles, with PERMANOVA, ANOSIM and PERMDISP.

Joint PCoA of the true profiles of held-out genomes and each method's predicted profiles for the same
genomes: KO presence (Jaccard) and module completeness (Bray-Curtis). Per method: Mantel r between
true and predicted genome-to-genome distances, Procrustes m2 of the 10-axis ordinations, and paired
group tests of predicted vs true profiles (labels swapped within each genome under the null).
"""
import numpy as np
import pandas as pd
from scipy.linalg import eigh
from scipy.spatial import procrustes
from scipy.spatial.distance import cdist
from scipy.stats import rankdata

from benchlib import load_pred


def jaccard_dist(A):
    A = A.astype("float32")
    inter = A @ A.T
    s = A.sum(1)
    union = s[:, None] + s[None, :] - inter
    return np.where(union > 0, 1 - inter / np.maximum(union, 1), 0.0).astype("float64")


def bray_curtis_dist(M):
    s = M.sum(1)
    return cdist(M, M, "cityblock") / np.maximum(s[:, None] + s[None, :], 1e-12)


def pcoa(D, k=2):
    """Classical scaling: top-k coordinates and their share of total inertia (trace of the centred matrix)."""
    D2 = D ** 2
    B = -0.5 * (D2 - D2.mean(0)[None, :] - D2.mean(1)[:, None] + D2.mean())
    n = len(B)
    w, v = eigh(B, subset_by_index=[n - k, n - 1])
    w, v = w[::-1], v[:, ::-1]
    return v * np.sqrt(np.maximum(w, 0)), w / np.trace(B)


def group_tests(D, n, rng, perms=999):
    """PERMANOVA, ANOSIM and PERMDISP for two groups of n paired profiles (rows 0..n-1 vs n..2n-1, row i
    paired with n+i: same genome). Permutations swap the two labels within each genome (strata = genome).
    Returns effect sizes, permutation p-values and each group's mean distance to its centroid."""
    N = 2 * n
    Q = D ** 2
    r = rankdata(D[np.triu_indices(N, 1)])
    R = np.zeros_like(D)
    R[np.triu_indices(N, 1)] = r
    R = R + R.T
    swap = rng.random((perms, n)) < 0.5
    A = np.zeros((N, perms + 1))
    A[n:, 0] = 1  # observed: second block is group 1
    A[:n, 1:] = swap.T
    A[n:, 1:] = ~swap.T
    n1 = n0 = n
    ones_q, ones_r = Q.sum(1), R.sum(1)

    def within(M, rowsum):
        MA = M @ A
        aMa = (A * MA).sum(0)
        bMb = rowsum.sum() - 2 * (rowsum @ A) + aMa
        return MA, aMa, bMb
    QA, aQa, bQb = within(Q, ones_q)
    ss_t = Q.sum() / (2 * N)
    ss_w = aQa / (2 * n1) + bQb / (2 * n0)
    F = (ss_t - ss_w) / (ss_w / (N - 2))
    _, aRa, bRb = within(R, ones_r)
    w_pairs = n1 * (n1 - 1) / 2 + n0 * (n0 - 1) / 2
    rw = (aRa + bRb) / 2 / w_pairs
    rb = (r.sum() - (aRa + bRb) / 2) / (n1 * n0)
    Ran = (rb - rw) / (len(r) / 2)
    QB = ones_q[:, None] - QA  # Q @ (1 - a)
    z2 = A * (QA / n1 - aQa / (2 * n1 ** 2)) + (1 - A) * (QB / n0 - bQb / (2 * n0 ** 2))
    z = np.sqrt(np.maximum(z2, 0))
    m1, m0 = (A * z).sum(0) / n1, ((1 - A) * z).sum(0) / n0
    zbar = z.mean(0)
    ssb = n1 * (m1 - zbar) ** 2 + n0 * (m0 - zbar) ** 2
    ssw = (A * (z - m1) ** 2).sum(0) + ((1 - A) * (z - m0) ** 2).sum(0)
    Fd = ssb / (ssw / (N - 2))
    # ties count as >= (tolerance: identical groups give equal statistics up to rounding)
    pval = lambda x: (1 + (x[1:] >= x[0] - 1e-9 * max(1.0, abs(x[0]))).sum()) / (perms + 1)  # noqa: E731
    return {"permanova_R2": (ss_t - ss_w[0]) / ss_t, "permanova_p": pval(F), "anosim_R": Ran[0],
            "anosim_p": pval(Ran), "dispersion_group0": m0[0], "dispersion_group1": m1[0],
            "permdisp_p": pval(Fd)}


def collect(tasks, methods, level, arm="full"):
    """Calls of each method for every held-out genome of one level x arm. methods: name -> (file stem,
    threshold or None: None = the file's own calls, a number = score >= threshold)."""
    calls, queries = {m: [] for m in methods}, []
    for t in tasks:
        info = pd.read_csv(t / "task.tsv", sep="\t")
        if info.empty or info.level.iloc[0] != level or info.arm.iloc[0] != arm:
            continue
        if not all((t / f"pred_{stem}.npz").exists() for stem, _ in methods.values()):
            continue
        queries += info.accession.tolist()
        for m, (stem, th) in methods.items():
            p = load_pred(t / f"pred_{stem}.npz")
            calls[m].append(p["call"] if th is None else p["score"] >= th)
    if not queries:
        return None, None
    return queries, {m: np.vstack(v) for m, v in calls.items()}


def ordination(tasks, truth_Y, truth_mod, scorer, rng, max_n, level, methods):
    """Joint PCoA per space, agreement with truth and paired group tests, for a subsample of genomes."""
    queries, calls = collect(tasks, methods, level)
    if queries is None:
        return None
    idx = np.arange(len(queries))
    if len(idx) > max_n:  # shortcut: dense n^2 PCoA and permutation tests; a random subsample of genomes
        idx = np.sort(rng.choice(idx, max_n, replace=False))
    q = [queries[i] for i in idx]
    kos = truth_Y.columns
    spaces = []
    for label, dist, truth, pred in [
        ("KO presence (Jaccard)", jaccard_dist, truth_Y.loc[q].to_numpy(),
         {m: c[idx] for m, c in calls.items()}),
        ("Module completeness (Bray-Curtis)", bray_curtis_dist, truth_mod.loc[q].to_numpy(),
         {m: scorer.matrix(pd.DataFrame(c[idx], index=q, columns=kos)).reindex(
             columns=truth_mod.columns, fill_value=0).to_numpy() for m, c in calls.items()})]:
        names = ["truth"] + list(pred)
        D = dist(np.vstack([truth] + [pred[m] for m in pred]))
        coords, expl = pcoa(D, 2)
        n = len(q)
        Dt = D[:n, :n]
        iu = np.triu_indices(n, 1)
        t10, _ = pcoa(Dt, 10)
        stats = []
        for k, m in enumerate(names[1:], start=1):
            blk = slice(k * n, (k + 1) * n)
            Dm = D[blk, blk]
            both = np.r_[np.arange(n), np.arange(k * n, (k + 1) * n)]
            gt = group_tests(D[np.ix_(both, both)], n, rng)
            stats.append({"space": label, "method": m, "genomes": n,
                          "mantel_r": np.corrcoef(Dt[iu], Dm[iu])[0, 1],
                          "procrustes_m2": procrustes(t10, pcoa(Dm, 10)[0])[2],
                          "distance_to_own_truth": np.diag(D[blk, :n]).mean(),
                          "permanova_R2": gt["permanova_R2"], "permanova_p": gt["permanova_p"],
                          "anosim_R": gt["anosim_R"], "anosim_p": gt["anosim_p"],
                          "spread_truth": gt["dispersion_group0"], "spread_method": gt["dispersion_group1"],
                          "permdisp_p": gt["permdisp_p"]})
        seg = rng.choice(n, min(120, n), replace=False)
        spaces.append(dict(label=label, coords=coords, expl=expl, n=n, names=names, stats=pd.DataFrame(stats), seg=seg))
    return dict(spaces=spaces, stats=pd.concat([s["stats"] for s in spaces], ignore_index=True), level=level)

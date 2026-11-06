"""Simulated microbial communities with exact functional truth, shared by the benchmark and the exploration.

Design (DESIGN §B.5, D095). A community is built only from genomes held out of every reference, one ASV per
genome (its 16S or amplicon). Each replicate draws a habitat pool; each sample holds a random subset of the pool
with log-normal cell abundances; in group B a set of responder genomes is `fold_change` times more abundant.
Reads per genome = cells x true 16S copies. A method turns reads back into cells with its predicted 16S copy
number and into KO abundance with its predicted KO copies, as PICRUSt2's metagenome step does.

Habitat designs approximate published V4 surveys (soil > marine > gut in richness and evenness; Earth
Microbiome Project, Shaffer et al. 2020, PMC7510982): richness = ASVs per sample, evenness = Pielou's J,
reached by calibrating the log-normal spread. These are design parameters, not measurements; the report
checks that conclusions hold across all three.
"""
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

HABITATS = {"gut": dict(richness=200, evenness=0.65),
            "marine": dict(richness=500, evenness=0.75),
            "soil": dict(richness=1000, evenness=0.88)}


def pielou(x):
    p = x[x > 0] / x[x > 0].sum()
    return -(p * np.log(p)).sum() / np.log(len(p)) if len(p) > 1 else 1.0


def calibrate_sigma(richness, evenness, rng, n=200):
    """Log-normal sigma whose expected Pielou J (for this richness) matches `evenness` (bisection)."""
    lo, hi = 0.05, 5.0
    for _ in range(30):
        mid = (lo + hi) / 2
        j = np.mean([pielou(rng.lognormal(0, mid, richness)) for _ in range(200)])
        lo, hi = (mid, hi) if j > evenness else (lo, mid)
    return (lo + hi) / 2


def simulate(n_genomes, habitat, rng, samples=100, responder_share=0.1, fold_change=4.0, occupancy=0.5):
    """Cells (samples x n_genomes, zeros for absent genomes), group (bool: B), responders (bool per genome)."""
    h = HABITATS[habitat]
    pool_n = min(n_genomes, int(round(h["richness"] / occupancy)))
    sigma = calibrate_sigma(h["richness"], h["evenness"], np.random.default_rng(0))
    pool = rng.choice(n_genomes, pool_n, replace=False)
    responders = np.zeros(n_genomes, bool)
    responders[rng.choice(pool, max(1, int(round(responder_share * pool_n))), replace=False)] = True
    group = np.arange(samples) >= samples // 2
    cells = np.zeros((samples, n_genomes))
    for s in range(samples):
        present = pool[rng.random(pool_n) < occupancy]
        cells[s, present] = rng.lognormal(0, sigma, len(present))
    cells[np.ix_(group, responders)] *= fold_change
    return cells, group, responders


def profile(cells, copies):
    """Relative function abundance per sample: cells @ copies, normalised to sum 1."""
    X = cells @ copies
    return X / np.maximum(X.sum(1, keepdims=True), 1e-12)


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    q = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1)
    return out


def da_calls(X, group):
    """Samples x functions -> BH-significant functions (Wilcoxon rank-sum A vs B, q < 0.05)."""
    p = mannwhitneyu(X[group], X[~group], axis=0).pvalue  # vectorised; same p-values as one test per column
    return bh(np.where(X.std(0) > 0, np.nan_to_num(p, nan=1.0), 1.0)) < 0.05


def log2_fold_change(X, group, eps=1e-9):
    return np.log2(X[group].mean(0) + eps) - np.log2(X[~group].mean(0) + eps)


def score(pred, truth, group):
    """Community metrics of one predicted profile table against the truth (samples x functions)."""
    bc = np.abs(pred - truth).sum(1) / np.maximum((pred + truth).sum(1), 1e-12)
    sp = [spearmanr(pred[i], truth[i]).statistic for i in range(len(truth))]
    keep = truth.std(0) > 0
    t_da, p_da = da_calls(truth[:, keep], group), da_calls(pred[:, keep], group)
    tp, fp, fn = (p_da & t_da).sum(), (p_da & ~t_da).sum(), (~p_da & t_da).sum()
    return dict(bray_curtis=float(np.median(bc)), spearman=float(np.median(sp)),
                lfc_spearman=float(spearmanr(log2_fold_change(pred[:, keep], group),
                                             log2_fold_change(truth[:, keep], group)).statistic),
                da_true=int(t_da.sum()), da_predicted=int(p_da.sum()), da_precision=tp / max(tp + fp, 1),
                da_recall=tp / max(tp + fn, 1), da_f1=2 * tp / max(2 * tp + fp + fn, 1))

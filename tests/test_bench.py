import numpy as np


def test_kmer_features_counts():
    from bench_mlp import kmer_features
    F = kmer_features(["ACGTA", "AANA"], kmax=2)
    assert F.shape == (2, 4 + 16)
    assert np.allclose(F[0, :4], [2 / 5, 1 / 5, 1 / 5, 1 / 5])
    # 2-mers of ACGTA: AC, CG, GT, TA; AANA: only AA (windows with N skipped)
    two = F[:, 4:]
    assert np.isclose(two[0].sum(), 1) and two[0, 0 * 4 + 1] == 0.25 and two[1, 0] == 1.0


def test_cluster_bootstrap_resamples_whole_clusters():
    from benchlib import cluster_bootstrap
    rng = np.random.default_rng(0)
    # two clusters: all +1 in one, all -1 in the other; resampling clusters gives a wide CI around 0
    d = np.r_[np.ones(500), -np.ones(500)]
    mean, lo, hi = cluster_bootstrap(d, ["a"] * 500 + ["b"] * 500, rng, 2000)
    assert mean == 0 and lo <= -0.99 and hi >= 0.99
    # independent genomes: narrow CI
    mean, lo, hi = cluster_bootstrap(d, np.arange(1000), rng, 2000)
    assert hi - lo < 0.2


def test_expected_f1_calls_picks_top_n():
    from s16fun import expected_f1_calls
    P = np.array([[0.95, 0.9, 0.6, 0.05, 0.01], [0.1, 0.1, 0.1, 0.1, 0.1]], np.float32)
    C = expected_f1_calls(P)
    assert C[0].tolist() == [True, True, True, False, False]
    assert C[1].sum() >= 1   # always calls at least one KO


def test_pred_roundtrip(tmp_path):
    from benchlib import load_pred, save_pred
    S = np.random.default_rng(1).random((3, 11)).astype(np.float32)
    save_pred(tmp_path / "p.npz", ["a", "b", "c"], S, S > 0.5, copies=S * 2)
    p = load_pred(tmp_path / "p.npz")
    assert p["queries"] == ["a", "b", "c"] and (p["call"] == (S > 0.5)).all()
    assert np.allclose(p["score"], S, atol=1e-3) and np.allclose(p["copies"], 2 * S, atol=2e-3)


def test_picrust2_ref_subset_excludes_held_out(tmp_path):
    import gzip
    from bench_picrust2_ref import build
    psc = tmp_path / "psc"
    for dom, ref in [("bacteria", "bac_ref"), ("archaea", "arc_ref")]:
        d = psc / dom / ref
        d.mkdir(parents=True)
        ids = [f"{ref}_{i}" for i in range(4)]
        (d / f"{ref}.fna").write_text("".join(f">{i}\nAC-GT\n" for i in ids))
        (d / f"{ref}.tre").write_text(f"(({ids[0]}:1,{ids[1]}:1):1,({ids[2]}:1,{ids[3]}:1):1);")
        (d / f"{ref}.model").write_text("GTR")
        with gzip.open(psc / dom / "ko.txt.gz", "wt") as f:
            f.write("assembly\tko:K1\tko:K2\n" + "".join(f"{i}\t1\t{int(n == 3)}\n" for n, i in enumerate(ids)))
        with gzip.open(psc / dom / "16S.txt.gz", "wt") as f:
            f.write("assembly\t16S_rRNA_Count\n" + "".join(f"{i}\t1\n" for i in ids))
    keep = {"bac_ref_0", "bac_ref_1", "bac_ref_2", "arc_ref_0", "arc_ref_1"}
    build(psc, keep, tmp_path / "out")
    fna = (tmp_path / "out/bac_ref/bac_ref.fna").read_text()
    assert "bac_ref_3" not in fna and "bac_ref_3" not in (tmp_path / "out/bac_ref/bac_ref.tre").read_text()
    with gzip.open(tmp_path / "out/bac_ko.tsv.gz", "rt") as f:
        lines = f.read().splitlines()
    assert lines[0] == "assembly" and len(lines) == 4   # K1 constant (1), K2 constant (0) once the held-out genome is gone
    const = (tmp_path / "out/bac_constant.tsv").read_text().splitlines()
    assert const == ["ko\tvalue", "ko:K1\t1", "ko:K2\t0"]


def test_s16fun_predicts_own_reference_exactly(tmp_path):
    """A query identical to a reference genome's 16S gets that genome's KOs (sanity check)."""
    import shutil
    import pandas as pd
    from scipy import sparse
    import s16fun as sf
    if shutil.which("vsearch") is None:
        import pytest
        pytest.skip("vsearch not installed")
    rng = np.random.default_rng(3)
    seqs = {f"g{i}": "".join(rng.choice(list("ACGT"), 600)) for i in range(30)}
    g = pd.DataFrame({r: ["x"] * 30 for r in sf.PRIOR_RANKS} | {"ssu_copies": [1] * 30}, index=sorted(seqs))
    g["genus"] = [f"G{i // 3}" for i in range(30)]
    Y = (rng.random((30, 40)) < 0.4).astype(np.int16)
    (tmp_path / "ref.fna").write_text("".join(f">{k}\n{seqs[k]}\n" for k in g.index))
    ref = sf.Reference(g, tmp_path / "ref.fna", [f"K{i}" for i in range(40)], sparse.csr_matrix(Y))
    (tmp_path / "q.fna").write_text(f">q\n{seqs['g7']}\n")
    r = sf.predict(ref, tmp_path / "q.fna", 2, tmp_path)
    assert r["novelty"].loc["q", "nearest"] == "g7" and r["novelty"].loc["q", "top_id"] == 1.0
    assert (r["ev"]["knn1"][0] == Y[list(g.index).index("g7")]).all()


def test_weighted_auroc_matches_sklearn():
    from sklearn.metrics import roc_auc_score
    from benchlib import WeightedAUROC
    rng = np.random.default_rng(5)
    Y = rng.random((60, 7)) < 0.4
    S = np.round(rng.random((60, 7)) * 3) / 3 + Y * 0.3      # ties on purpose
    w = rng.integers(0, 4, 60).astype(float)
    got = WeightedAUROC(S, Y)(w)
    for j in range(7):
        assert np.isclose(got[j], roc_auc_score(Y[:, j], S[:, j], sample_weight=w))
    assert np.allclose(WeightedAUROC(S, Y)(np.ones(60)), [roc_auc_score(Y[:, j], S[:, j]) for j in range(7)])


def test_weighted_median_of_neighbours():
    from s16fun import weighted_median
    vals = np.array([1.0, 2.0, 7.0, 4.0])
    nbr = np.array([[0, 1, 2, -1, -1], [3, -1, -1, -1, -1], [-1, -1, -1, -1, -1]])
    w = np.array([[1, 1, 5, 0, 0], [1, 0, 0, 0, 0], [0, 0, 0, 0, 0]], float)
    assert weighted_median(vals, nbr, w, 5).tolist() == [7.0, 4.0, 1.0]

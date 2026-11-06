import numpy as np
from scipy.spatial.distance import pdist, squareform

from picrust2_compare import bray_curtis_dist, jaccard_dist, ko_scores, pcoa
from picrust2_eval import in_reference


def test_distances_match_scipy():
    rng = np.random.default_rng(0)
    A = rng.random((30, 50)) < 0.3
    assert np.allclose(jaccard_dist(A), squareform(pdist(A, "jaccard")))
    M = rng.random((30, 12))
    assert np.allclose(bray_curtis_dist(M), squareform(pdist(M, "braycurtis")))


def test_pcoa_recovers_planar_configuration():
    X = np.random.default_rng(1).normal(size=(40, 2)) * [3, 1]
    D = squareform(pdist(X))
    coords, expl = pcoa(D, 2)
    assert np.allclose(squareform(pdist(coords)), D, atol=1e-8)
    assert np.isclose(expl.sum(), 1) and expl[0] > expl[1]


def test_ko_scores_and_reference_flag(tmp_path):
    Y = np.array([[1, 1, 0, 0]], bool)
    P = np.array([[1, 0, 1, 0]], bool)
    s = ko_scores(Y, P)
    assert s["ko_f1"][0] == 0.5 and s["ko_precision"][0] == 0.5 and s["ko_recall"][0] == 0.5
    refs = tmp_path / "refs.txt"
    refs.write_text("GCA_000005845.1\nGCF_000001.2\n")
    flags = in_reference(["RS_GCF_000005845.2", "GB_GCA_000009999.1", "RS_GCF_000001.1"], refs)
    assert flags.tolist() == [True, False, True]


def test_group_tests_null_and_effect():
    from picrust2_compare import group_tests
    rng = np.random.default_rng(0)
    X = rng.random((150, 40)) < 0.3
    same = group_tests(jaccard_dist(np.vstack([X, X])), 150, rng)
    assert same["permanova_p"] == 1 and same["permdisp_p"] == 1 and same["anosim_p"] == 1
    null = group_tests(jaccard_dist(np.vstack([X, rng.random((150, 40)) < 0.3])), 150, rng)
    assert null["permanova_p"] > 0.01 and null["permanova_R2"] < 0.01
    eff = group_tests(jaccard_dist(np.vstack([X, rng.random((150, 40)) < 0.6])), 150, rng)
    assert eff["permanova_p"] <= 0.002 and eff["anosim_R"] > 0.1 and eff["dispersion_group1"] < eff["dispersion_group0"]

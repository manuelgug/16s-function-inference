import numpy as np
import pandas as pd
import pytest

from aim1_predict import base_onehot, genome_metrics, identity, knn_predict


def test_identity_ignores_gaps_and_counts_matches():
    chars = pd.DataFrame([list("ACGT"), list("ACGA"), list("AC--")])
    oh = base_onehot(chars)
    sim = identity(oh, oh)
    assert sim[0, 0] == 1 and sim[0, 1] == pytest.approx(0.75)
    assert sim[0, 2] == pytest.approx(0.5)  # 2 matches / 4 bases in the longer sequence


def test_knn1_copies_nearest_and_knn5_weights():
    sim = np.array([[0.9, 0.5]])
    Y = np.array([[1.0, 0.0], [0.0, 1.0]])
    assert knn_predict(sim, Y, 1).tolist() == [[1.0, 0.0]]
    p = knn_predict(sim, Y, 2)
    assert p[0, 0] > p[0, 1] and p.sum() == pytest.approx(1.0)


def test_ko_f1_and_jaccard():
    Y = np.array([[1, 1, 0, 0]], dtype=float)
    P = np.array([[0.9, 0.2, 0.7, 0.1]])
    m = genome_metrics("ko", Y, P)
    assert m["f1"][0] == pytest.approx(0.5) and m["jaccard"][0] == pytest.approx(1 / 3)


def test_gain_sign_lower_is_better():
    from aim1_report import gains
    g = pd.DataFrame({"arm": "full", "level": "random", "accession": "a", "target_type": "module", "metric": "mae",
                      "model": ["M0_null", "M3_rf"], "value": [0.20, 0.15], "fold": 0})
    assert gains(g)["gain"].iloc[0] == pytest.approx(0.05)  # lower MAE than null = positive gain

from pathlib import Path

import pytest

from build_targets import cog_matrix, ko_matrix, ko_sets, read_annotations

ANN = """## emapper-2.1.15
## command: emapper.py ...
#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko
gA@1\tx\t0\t1\t-\t-\tC\t-\t-\t-\t-\tko:K00844,ko:K01810
gA@2\tx\t0\t1\t-\t-\tEG\t-\t-\t-\t-\t-
gB@1\tx\t0\t1\t-\t-\t-\t-\t-\t-\t-\tko:K00844
## 3 queries scanned
"""


@pytest.fixture
def ann(tmp_path):
    p = tmp_path / "x.emapper.annotations"
    p.write_text(ANN)
    return read_annotations([p])


def test_read_and_ko_matrix(ann):
    m = ko_matrix(ko_sets(ann), ["gA", "gB"])
    assert list(m.columns) == ["K00844", "K01810"]
    assert m.loc["gA"].tolist() == [1, 1] and m.loc["gB"].tolist() == [1, 0]


def test_cog_fractions_split_multi_letter(ann):
    m = cog_matrix(ann, ["gA", "gB"])
    assert m.loc["gA", "C"] == pytest.approx(0.5)
    assert m.loc["gA", "E"] == pytest.approx(0.25) and m.loc["gA", "G"] == pytest.approx(0.25)
    assert m.loc["gB"].sum() == 0  # no COG-annotated proteins

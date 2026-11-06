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


def test_module_scorer_matches_give_completeness(tmp_path):
    """Vectorised completeness must equal the library's command-line tool, genome by genome."""
    import subprocess

    import numpy as np
    import pandas as pd

    from build_targets import ModuleScorer
    ms = ModuleScorer()
    rng = np.random.default_rng(0)
    pool = sorted({k for labels, _, _ in ms.paths.values() for k in labels if k.startswith("K")})
    genomes = {f"g{i}": set(rng.choice(pool, size=rng.integers(1, 400), replace=False)) for i in range(6)}
    genomes["g_none"] = {pool[0]}
    inp = tmp_path / "in.tsv"
    inp.write_text("".join(g + "\t" + "\t".join(sorted(k)) + "\n" for g, k in genomes.items()))
    subprocess.run(["give_completeness", "-i", str(inp), "-m", "-o", str(tmp_path / "o")], check=True,
                   stdout=subprocess.DEVNULL)
    ref = pd.read_csv(tmp_path / "o" / "summary.kegg_contigs.tsv", sep="\t")
    ref = ref.pivot(index="contig", columns="module_accession", values="completeness").fillna(0) / 100
    presence = pd.DataFrame([{k: 1 for k in v} for v in genomes.values()], index=list(genomes)).fillna(0)
    new = ms.matrix(presence)
    ref = ref.reindex(index=new.index, columns=new.columns).fillna(0)
    assert set(ref.columns[ref.any()]) <= set(new.columns)
    assert np.allclose(new.to_numpy(), ref.to_numpy(), atol=1e-6)

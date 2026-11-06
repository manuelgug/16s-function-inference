import numpy as np
import pandas as pd

from make_splits import RANKS, make_splits


def toy(n=400, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        p, c, o, f, g = (rng.integers(3), rng.integers(3), rng.integers(3), rng.integers(4), rng.integers(4))
        rows.append(dict(accession=f"G{i:04d}", domain="Bacteria", phylum=f"p{p}", **{"class": f"c{p}{c}"},
                         order=f"o{p}{c}{o}", family=f"f{p}{c}{o}{f}", genus=f"g{p}{c}{o}{f}{g}",
                         species=f"s{i}", group=f"G{i:04d}"))
    df = pd.DataFrame(rows).set_index("accession")
    df.loc["G0001", "group"] = "G0000"   # identical 16S pair across genera
    return df


def test_every_genome_tested_once_per_level():
    s = make_splits(toy(), ["random", "genus", "family"], 5, 1)
    assert (s.groupby("level").accession.nunique() == 400).all()
    assert (s.groupby(["level", "accession"]).size() == 1).all()


def test_holdout_never_splits_a_taxon_or_a_16s_group():
    df = toy()
    s = make_splits(df, ["random", "genus", "family", "order", "class"], 5, 1)
    for level in ["genus", "family", "order", "class"]:
        sl = s[s.level == level].set_index("accession")
        first = df.groupby("group").head(1).set_index("group")[level]
        taxon = first.loc[df.group].set_axis(df.index)
        assert (sl.groupby(taxon.loc[sl.index]).fold.nunique() == 1).all()
    for level in s.level.unique():
        sl = s[s.level == level].set_index("accession")
        assert sl.loc["G0000", "fold"] == sl.loc["G0001", "fold"]


def test_partition_is_by_family_and_same_for_every_level():
    df = toy()
    s = make_splits(df, ["random", "genus", "class"], 5, 1)
    part = s.pivot(index="accession", columns="level", values="partition")
    assert (part.nunique(axis=1) == 1).all()
    p = part["genus"]
    assert (p.groupby(df.loc[p.index, "family"]).nunique() == 1).all()
    assert 0.3 < (p == "dev").mean() < 0.7


def test_deterministic():
    assert make_splits(toy(), ["family"], 5, 3).equals(make_splits(toy(), ["family"], 5, 3))

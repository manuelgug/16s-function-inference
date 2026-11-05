import pandas as pd

from make_splits import make_splits


def toy():
    rows = []
    for fam in range(10):
        for gen in range(3):
            for sp in range(4):
                rows.append({"accession": f"f{fam}g{gen}s{sp}", "family": f"F{fam}", "genus": f"F{fam}G{gen}"})
    return pd.DataFrame(rows).set_index("accession")


def test_every_genome_tested_once_per_level():
    s = make_splits(toy(), ["random", "genus", "family"], 5, seed=1)
    counts = s.groupby(["level", "accession"]).size()
    assert (counts == 1).all() and len(counts) == 3 * 120


def test_group_levels_never_split_a_taxon():
    df = toy()
    s = make_splits(df, ["genus", "family"], 5, seed=1).join(df, on="accession")
    for level in ["genus", "family"]:
        assert (s[s.level == level].groupby(level)["fold"].nunique() == 1).all()


def test_deterministic():
    a = make_splits(toy(), ["random", "family"], 5, seed=7)
    b = make_splits(toy(), ["random", "family"], 5, seed=7)
    assert a.equals(b)

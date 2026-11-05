import pandas as pd
import pytest

from select_genomes import allocate, stratified_sample, subsample


def toy(sizes):
    rows = []
    for p, k in sizes.items():
        for i in range(k):
            rows.append({"accession": f"{p}_{i:04d}", "phylum": p,
                         "domain": "Archaea" if p.startswith("A") else "Bacteria"})
    return pd.DataFrame(rows)


def test_allocate_respects_floor_total_and_caps():
    sizes = {"big": 1000, "mid": 100, "tiny": 3}
    a = allocate(sizes, 200, floor=10)
    assert a.sum() == 200
    assert a["tiny"] == 3  # capped at availability
    assert a["mid"] >= 10 and a["big"] > a["mid"]


def test_allocate_rejects_n_below_floors():
    with pytest.raises(ValueError):
        allocate({"a": 50, "b": 50}, 5, floor=10)


def test_stratified_sample_deterministic_and_sized():
    df = toy({"big": 500, "mid": 60, "tiny": 2})
    s1 = stratified_sample(df, 100, 10, seed=1)
    s2 = stratified_sample(df, 100, 10, seed=1)
    assert len(s1) == 100 and s1.accession.is_unique
    assert list(s1.accession) == list(s2.accession)
    assert (s1.phylum == "tiny").sum() == 2


def test_stratified_sample_n0_keeps_all():
    df = toy({"x": 5})
    assert len(stratified_sample(df, 0, 10, seed=1)) == 5


def test_subsample_guarantees_archaea():
    df = toy({"Arc": 5, "Bac": 100})
    s = subsample(df, 10, min_archaea=2, seed=3)
    assert len(s) == 10 and (s.domain == "Archaea").sum() >= 2

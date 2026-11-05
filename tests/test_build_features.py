import numpy as np
import pandas as pd

from build_features import allele_features, ecoli_positions, main, read_alignment

# 4 genomes + reference; column 3 is a gap in the reference (insertion), column 5 has a SNP
ALN = """>ecoli_rrsB
GTGC-CAGCAAAATTTTCCCCGGGGATTAGATACCC
>g1
GTGCACAGCAAAATTTTCCCCGGGGATTAGATACCC
>g2
GTGCACAGCAAAATTTTCCCCGGGGATTAGATACCC
>g3
GTGCATAGCAAAATTTTCCCCGGGGATTAGATACCC
>g4
GTGCATAGCAAAATTTTCCCCGGGGATTAGATACCC
"""


def test_ecoli_positions_and_alleles(tmp_path):
    p = tmp_path / "a.fasta"
    p.write_text(ALN)
    ids, aln = read_alignment(p)
    pos = ecoli_positions(aln[0])
    assert pos[3] == 4 and pos[4] == 4 and pos[5] == 5  # insertion shares E. coli number
    df, meta = allele_features(aln[1:], ids[1:], np.arange(aln.shape[1]), pos, 0.01)
    assert set(df.columns) == {"5_5_C", "5_5_T"}
    assert df["5_5_T"].tolist() == [0, 0, 1, 1]


def test_main_with_pcr_arm(tmp_path):
    (tmp_path / "a.fasta").write_text(ALN)
    (tmp_path / "ssu.fna").write_text("".join(
        f">{g}\n{s.replace('-', '')}\n" for g, s in
        zip(["g1", "g2", "g3", "g4"], [l for l in ALN.splitlines()[3::2]])))
    # forward primer matches the first 4 bases, reverse primer the last 6 (reverse complemented)
    (tmp_path / "primers.tsv").write_text("toy\tGTGC\tGGGTAT\n")
    main(["--aln", str(tmp_path / "a.fasta"), "--ssu", str(tmp_path / "ssu.fna"),
          "--primers", str(tmp_path / "primers.tsv"), "--outdir", str(tmp_path / "o")])
    full = pd.read_parquet(tmp_path / "o" / "features_full.parquet")
    arm = pd.read_parquet(tmp_path / "o" / "features_toy.parquet")
    pcr = pd.read_csv(tmp_path / "o" / "pcr_toy.tsv", sep="\t")
    assert pcr.pcr.sum() == 4
    assert list(arm.columns) == list(full.columns)  # SNP lies inside the amplicon

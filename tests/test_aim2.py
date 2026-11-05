import dendropy
import numpy as np

from aim2_prep import kinship


def test_kinship_is_shared_root_to_mrca_length():
    tree = dendropy.Tree.get(data="((A:1,B:2):3,C:4);", schema="newick", rooting="force-rooted")
    K = kinship(tree, ["A", "B", "C"])
    assert np.allclose(K, [[4, 3, 0], [3, 5, 0], [0, 0, 4]])


def test_fitch_changes():
    import pandas as pd

    from aim2_report import fitch_changes
    tree = dendropy.Tree.get(data="((A:1,B:1):1,(C:1,D:1):1);", schema="newick", rooting="force-rooted")
    geno = pd.DataFrame({"homoplasic": [1, 0, 1, 0], "clade": [1, 1, 0, 0], "const": [0, 0, 0, 0]},
                        index=list("ABCD"))
    assert fitch_changes(tree, geno).to_dict() == {"homoplasic": 2, "clade": 1, "const": 0}

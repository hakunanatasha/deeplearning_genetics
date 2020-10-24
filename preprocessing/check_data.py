"""
Natasha Seelam (nseelam1@gmail.com)
2020.10.24

The following checks whether your dataset
exhibits the A/B motifs of interest.

This is specific to a 0/1/2 motif situation, where
we look at data where no A or no B, only A, only B, or both A and B are present
in the randomly generated DNA sequences.

This can serve as a template for more complicated motif creation.
"""

import numpy as np
import pandas as pd

motifA = "CGACCGAACTCC"
motifB = "ACATGCTTAGTA"

savename_noAB = "../data/motif_noA_noB.csv"
savename_A = "../data/motif_A.csv"
savename_B = "../data/motif_B.csv"
savename_AB = "../data/motif_AB.csv"


# ---------------- #
# Load Dataset
# ---------------- #
noab = pd.read_csv(savename_noAB)
a = pd.read_csv(savename_A)
b = pd.read_csv(savename_B)
ab = pd.read_csv(savename_AB)


def check_counts(seqs, A, B):
    """
    Check if the number of motifs are present in the sequence
    """
    N = seqs.shape[0]
    t1 = [seqs.iloc[idx, 0].count(A) - seqs.iloc[idx, 1] for idx in range(N)]
    t2 = [seqs.iloc[idx, 0].count(B) - seqs.iloc[idx, 2] for idx in range(N)]
    return sum(t1), sum(t2)


for key, value in {"No A no B": noab, "Aonly": a, "Bonly": b, "ABonly": ab}.items():
    print("\n" + key)
    isA, isB = check_counts(value, motifA, motifB)
    print("Matched A:", isA == 0)
    print("Matched B:", isB == 0)

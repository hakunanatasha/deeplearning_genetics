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

def check_counts(seqs, motif, col_idx):
    """
    Check if the number of motifs are present in the sequence

    Args:
    seqs: [pd.DataFrame] pandas data frame of [seqs, number_motifs_i]
    motif: [str] the motif substr you want to string match
    col_idx: [int] the column index that counts the motifs
    """
    N = seqs.shape[0]
    t1 = [seqs.iloc[idx, col_idx].count(motif) - seqs.iloc[idx][col_idx] for idx in range(N)]
    return sum(t1)


if __name__ == "__main__":

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

    for key, value in {"No A no B": noab, "Aonly": a, "Bonly": b, "ABonly": ab}.items():
        print("\n" + key)
        isA = check_counts(value, motifA, 1)
        isB = check_counts(value, motifB, 2)
        print("Matched A:", isA == 0)
        print("Matched B:", isB == 0)

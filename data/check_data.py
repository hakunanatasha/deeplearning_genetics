import numpy as np
import pandas as pd

motifA = "CGACCGAACTCC"
motifB = "ACATGCTTAGTA"

noab = pd.read_csv("motif_noA_noB.csv")
a = pd.read_csv("motif_A.csv")
b = pd.read_csv("motif_B.csv")
ab = pd.read_csv("motif_AB.csv")


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

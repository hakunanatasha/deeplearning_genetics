"""
Data generator

The following creates hypothetical sequences of DNA with a
guaranteed motif of interest.

The user can provide the following to customize their sequences:

(A) Motif(s) sequences
(B) Alphabet (Standard ATCG or the IUPAC naming scheme)
(C) Frequency of letter appearance.
"""

import numpy as np
from random import shuffle
from itertools import zip_longest
import sys


# Generate with no motifs, generate with 1 motif A, 1 motif B, and both motifs


def generate_random_seq(alphabet, N):
    letters, probs = zip(*list(alphabet.items()))
    return "".join(np.random.choice(letters, p=probs, size=N))


def generate_data(N_seq, alphabet, motifs, morder, seqlen_min, seqlen_max):
    """
    Args:
    N_seq: [int] number of sequences to generate
    alphabet: [dict[str: int]] letters in sequence and probability of appearing.
    motifs: [dict[str: int]] motif_sequence: number_of_this_motif
    morder: [dict[int: str]] motif_index: motif_sequence
    seqlen_min: [int] minimum sequence length
    seqlen_max: [int] maximum sequence length
    """

    # Number of motifs
    Nmotifs = len(motifs)

    # Initialize sequences + labels
    seq_labels = np.zeros(shape=(N_seq, Nmotifs))
    seqs = []

    for idx in range(N_seq):
        # Length of generated sequence
        Ni_seqlen = np.random.randint(seqlen_min, seqlen_max, 1).item()

        # Generate the frequency of a motif appearing in the sequence
        # Generate a random number of non-zero motifs up to v
        mfreq = {
            k: np.random.randint(low=1, high=v) for k, v in motifs.items() if v > 1
        }

        # Generate all motifs that appear 0 or 1 times.
        mfreq.update({k: v for k, v in motifs.items() if k not in mfreq.keys()})

        # Create the exact sequence and count the number of times a motif appeared
        seqs.append(generate(alphabet, Ni_seqlen, mfreq))
        seq_labels[idx, :] = [mfreq[morder[k]] for k in range(Nmotifs)]

    return seqs, seq_labels


def generate(alphabet, N, motifs):
    """
    Generate a synthetic DNA with exactly the number of motifs specified

    Args:
    alphabet: [dict] letter: p_genletter
    N: [int] length of sequence
    motifs: dict: number of times each motif must appear (no more, no less)
    """

    # duplicate motifs by value
    flat_motifs = [_ for k, v in motifs.items() for _ in [k] * v]
    shuffle(flat_motifs)

    # A "segment" is a region of DNA that appears either between motifs or at
    # either end of the sequence (i.e. linker or end-cap).
    num_segments = len(flat_motifs) + 1

    num_bp_to_generate = N - sum(len(m) for m in flat_motifs)

    pvals = [1 / num_segments] * num_segments  # uniform probs
    segment_lens = np.random.multinomial(num_bp_to_generate, pvals)
    segments = [generate_random_seq(alphabet, L) for L in segment_lens]

    seq = "".join(s + m for s, m in zip_longest(segments, flat_motifs, fillvalue=""))

    seq_motifs = {m: seq.count(m) for m in motifs}
    if seq_motifs != motifs:
        # try again if sequence doesn't contain the precise number of motifs
        return generate(alphabet, N, motifs)
    else:
        return seq

    seqs = [[k for _ in range(v)] for k, v in motifs.items()]

    return seq


if __name__ == "__main__":

    # User Parameters
    # Max number of A/B motifs
    maxA = int(sys.argv[1])
    maxB = int(sys.argv[2])

    # Data CSV name
    savename = str(sys.argv[3])

    # Define letters and alphabet probability
    letters = ["A", "T", "C", "G"]
    alphabet = {l: 1 / len(letters) for l in letters}

    # Given motifs
    motifA = "CGACCGAACTCC"
    motifB = "ACATGCTTAGTA"

    # TODO: Later generalize to multiple motifs
    # with open('motif.txt', 'r') as f:
    #    motifA, motifB = [x.split()[0] for x in f.readlines()]

    Nseqs = int(1e4)

    # Length of the sequence min/max
    seqlen_min = int(1e3)
    seqlen_max = int(2e3)

    # Motif Theme, number of times motif appears, order
    mdict = {motifA: maxA, motifB: maxB}
    morder = {0: motifA, 1: motifB}

    print("Making", Nseqs, "Sequences with Motif A=", maxA, "|| Motif B=", maxB)

    motifs, labels = generate_data(
        Nseqs, alphabet, mdict, morder, seqlen_min, seqlen_max
    )

    with open(savename + ".csv", "w") as f:
        lines = ["Seq,A,B\n"]
        lines += [
            motifs[k] + "," + ", ".join([str(m) for m in labels[k, :]]) + "\n"
            for k in range(Nseqs)
        ]
        f.writelines(lines)

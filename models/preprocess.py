"""
Returns a dataset of one-hot-encoded sequences of DNA and whether the motif of interest is present or absent in the data.
"""

import pandas as pd
import numpy as np

from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.model_selection import train_test_split
import joblib as jb


class MakeOHE:
    """
    Make binary one-hot-encoders

    Args:
    motif_files: dict(label_id: list_of_data_filepaths)
    alphabet: list[str] Exhaustive list of all letters available in sequences
    """

    def __init__(self, motif_files, alphabet=["A", "C", "G", "T"]):
        self.motif_files = motif_files
        self.alphabet = alphabet

        self._load_data()
        self._get_length()
        self._setup_iencoder()
        self._setup_oheencoder()
        self.OHEncode()
        self.shuffle_data()

    def _load_data(self):
        print("Loading CSV")
        self.data = {
            k: pd.concat([pd.read_csv(fn) for fn in v])
            for k, v in self.motif_files.items()
        }
        self.label_ids = list(self.data.keys())

    def _get_length(self):
        self.N = {k: v.shape[0] for k, v in self.data.items()}

    def _setup_iencoder(self):
        print("Setting up Integer Encoder")
        self.i_encoder = LabelEncoder().fit(self.alphabet)

    def _setup_oheencoder(self, categories="auto"):
        print("Setting up OHE Encoder")
        self.input_encoder = OneHotEncoder(categories=categories)
        self.output_encoder = OneHotEncoder(categories=categories)

    def OHEncode(self):
        """
        Encode each data
        """
        print("Encoding Data.")
        data = []
        labels = []
        countr = 0
        for k, df in self.data.items():
            print("Encoded Classes =", countr + 1, "/", len(self.data))
            countr += 1
            # Get training sequences
            data += [
                self.ohe(df.iloc[idx, :], self.i_encoder, self.input_encoder)[0]
                for idx in range(self.N[k])
            ]
            # Get Labels
            labels += [k for _ in range(self.N[k])]

        self.xdata = data

        self.ydata = self.output_encoder.fit_transform(np.array(labels).reshape(-1, 1))
        print("Completed!")

    def shuffle_data(self, ptest=0.2, seed=1234):
        """
        Given a list of OHE training and OHE labels,
        shuffles the order.
        """
        print(
            "Creating shuffled training/testing data with ptest:",
            ptest,
            ", seed:",
            seed,
        )
        xtrain, xtest, ytrain, ytest = train_test_split(
            self.xdata, self.ydata, test_size=ptest, random_state=seed
        )

        self.xtrain = xtrain
        self.ytrain = ytrain
        self.xtest = xtest
        self.ytest = ytest

    @staticmethod
    def ohe(row, integer_encoder, one_hot_encoder):
        """
        One-hot-encodes A/T/C/G

        Returns:
        one_hot_encoded - tensor of Nseq x Nbases
        (intA, intB) - number of A and B
        """
        seq, intA, intB = row.tolist()

        # Find out which label position it is
        ienc = integer_encoder.transform(list(seq))
        ienc = np.array(ienc).reshape(-1, 1)

        # One-hot encode the position on an N x 4 vector
        one_hot_encoded = one_hot_encoder.fit_transform(ienc)
        return one_hot_encoded.toarray(), (intA, intB)

"""
Returns a dataset of one-hot-encoded sequences of DNA and whether the motif of interest is present or absent in the data.
"""

import pandas as pd
import numpy as np

# ML Preprocessing (OHE of data + model selection)
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.model_selection import StratifiedShuffleSplit

# Sequence padding
import torch
from torch.nn.utils.rnn import pad_sequence


class makeOHE:
    """
    Make binary one-hot-encoders

    Args:
    motif_files: dict(class_label: list_of_data_filepaths)
    alphabet: list[str] Exhaustive list of all letters available in sequences

    Takes the dictionary of class labels: list_of_CSV_files.
    Loads data, one-hot-encodes the (biological) sequences.
    Splits into train/test that are equally represented across classes.

    Args:
    motif_files: [dict class_label : path_to_class_data]
    alphabet: [list(str)] list of characters that appear in the sequences
    ptest: [float <1] percent of testing split
    pval: [float <1] percent of validation split
    batchfirst: [Bool] if data is returned as Nbatch x Nseq x Nfts
    process: [Bool] user has option to run the processing commands at init()
    """

    def __init__(
        self,
        motif_files,
        alphabet=["A", "C", "G", "T"],
        ptest=0.1,
        pval=0.2,
        seed=1234,
        batchfirst=True,
        pad_value=-1,
        process=True,
    ):
        assert (ptest + pval) < 1, "p_test + p_val are too big"

        self.motif_files = motif_files
        self.alphabet = sorted(alphabet)  # Alphabetically sort labels
        self.alph_dict = {a: idx for idx, a in enumerate(self.alphabet)}

        self.ptest = ptest
        self.pval = pval
        self.seed = seed
        self.batchfirst = batchfirst
        self.pad_value = -1

        if process:
            # Load all the CSVs and concatenate each class
            self.setup()

            # Encode each class and pad the sequences
            self.OHEncode()

            # Create a train/test split for each class
            print(
                "Train %=",
                (1 - self.ptest - self.pval),
                "|| Test %=",
                self.ptest,
                "|| Val %=",
                self.pval,
            )
            self.split_traintest()
        else:
            print("User did not specify processing steps.")

    def setup(self):
        """
        Loads motif data, gets the length and sets up one-hot-encodes
        """
        self._check_input_list()
        self._load_data()
        self._get_length()
        self._setup_encoders()

    def OHEncode(self):
        """
        Creates a one-hot-encoded vector for each sequence in all classes.
        - First, encodes all classes.
        - Then passes the full dataset to an RNN sequence padder that appends
        pad_value to all missing elements of the sequence and splits on class-size

        Returns: dict[Class_Label(str): torch.Tensor[Nbatch x Nseq_max x N_class]]
        """
        Nalphabet = len(self.alphabet)
        Nclasses = len(list(self.data.items()))

        data = {}
        labels = {}
        countr = 0
        for k, df in self.data.items():
            print("Encoding label=", k, "|| Classes =", countr + 1, "/", len(self.data))
            countr += 1
            # Get training sequences
            data.update(
                {
                    k: [
                        self.ohe(df.iloc[idx, :], self.i_encoder, Nalphabet)
                        for idx in range(self.N[k])
                    ]
                }
            )
            # Get Labels (list of the label repeating 'k' times)
            lval = [k for _ in range(self.N[k])]
            lval = self.ohe(lval, self.o_encoder, Nclasses)
            labels.update({k: lval})

        print("Padding")
        self.xdata = self.add_padding(data, self.N, self.batchfirst, self.pad_value)

        self.ydata = labels
        print("Completed!")

    def split_traintest(self):
        """
        Wrapper around stratified-shuffle-split.

        Ensures the training/testing split are equally represented
        in the input classes.

        Returns xtrain/test ytrain/test assigned
        where x represents the concatenated data and y the concatenated labels
        """
        print("Splitting Training + Testing")
        y = np.array([k for k, v in self.N.items() for _ in range(v)])
        X = torch.cat(list(self.xdata.values()))
        labels = torch.from_numpy(np.concatenate(list(self.ydata.values())))

        # First, create the Training/ [Validation + Testing] set
        stratSplit = StratifiedShuffleSplit(
            n_splits=1,
            test_size=round(self.ptest + self.pval, 2),
            random_state=self.seed,
        )
        for train_idx, out_idx in stratSplit.split(np.zeros(X.shape[0]), y):
            self.xtrain, xout = X[train_idx], X[out_idx]
            self.ytrain, yout = labels[train_idx], labels[out_idx]
            ytmp = y[out_idx]

        # From the reserved data, create validation/testing
        if self.pval != 0:
            print("Making Validation + Testing")
            ptest_2 = self.ptest / (self.ptest + self.pval)  # Scaled Ratio
            valSplit = StratifiedShuffleSplit(
                n_splits=1, test_size=ptest_2, random_state=self.seed
            )  # make a new split with remaining data
            for val_idx, test_idx in valSplit.split(np.zeros(xout.shape[0]), ytmp):
                self.xval, self.xtest = xout[val_idx], xout[test_idx]
                self.yval, self.ytest = yout[val_idx], yout[test_idx]

        else:
            self.xtest = xout
            self.ytest = yout
            self.xval = None
            self.yval = None

    def _check_input_list(self):
        """
        Input dictionary of motif files should be list. Converts.
        """
        f = lambda v: v if isinstance(v[1], list) else (v[0], [v[1]])
        self.motif_files = dict(map(f, self.motif_files.items()))

    def _load_data(self):
        """
        Load the sequence data from the CSV files of motif_files.
        Concatenates the set of CSV files.
        Expects all CSVs to have the format: Seq[str], Motif_1[int],...,Motif_i[int]
        where the column represents the number of times that motif appeared in the sequence.
        """
        print("Loading CSV")
        self.data = {
            k: pd.concat([pd.read_csv(fn) for fn in v])
            for k, v in self.motif_files.items()
        }
        self.label_ids = list(self.data.keys())

    def _get_length(self):
        """
        Get the number of data points.

        Returns dict[Class_Label(str): number_of_examples(int)]
        """
        self.N = {k: v.shape[0] for k, v in self.data.items()}

    def _setup_encoders(self):
        """
        Given the alphabet (ATCG or 20 amino acids) set up an encoder
        to account for each dimension.

        Then, set up the one-hot encoding s.t. when you see the label,
        you generate a "1" in the appropriate dimension.
        """
        print("Setting up Integer/Label Encoders")
        self.i_encoder = LabelEncoder().fit(self.alphabet)
        self.o_encoder = LabelEncoder().fit(list(self.data.keys()))

    @staticmethod
    def add_padding(data, N, batchfirst=True, pad_value=-1):
        """
        Pad all data sequences for each class.
        Collects all datapoints and then splits based on length of class
        """
        Ns = [0] + np.cumsum(list(N.values())).tolist()
        seqs = []
        pad = {}
        for k, v in data.items():
            seqs += v

        # Ensure all sequences are padded to same length
        alldata = pad_sequence([torch.Tensor(i) for i in seqs], batchfirst, pad_value)

        # Split on the classes
        Nsplits = [alldata[Ns[i] : Ns[i + 1], :, :] for i in range(len(N))]

        pad.update({k: Nsplits[idx] for idx, k in enumerate(data.keys())})

        return pad

    @staticmethod
    def ohe(row, integer_encoder, dimensions):
        """
        One-hot-encodes the alphabet.
        Row is a single row of a CSV file from self.data
        OR a list with the class of interest.

        Returns:
        one_hot_encoded - tensor of Nseq x Nbases
        (intA, intB) - number of A and B
        """
        if isinstance(row, list) is False:
            seq, intA, intB = row.tolist()
        else:
            seq = row
            intA = intB = None

        # Find out which label position it is
        ienc = integer_encoder.transform(list(seq))

        # One-hot encode the position on an N x 4 vector
        one_hot_encoded = np.zeros(shape=(len(seq), dimensions))
        one_hot_encoded[list(range(len(ienc))), ienc] = 1

        return one_hot_encoded  # , (intA, intB)

# deeplearning-genetics
An end-to-end model to learn genetic information for patient outcomes

## Data Generation [data/]

The following creates arbitrary genetic sequences with motif(s) of interest. User can specify the size of sequence, motif pattern, and the frequency of ATCG content. User can also specify what part of the IUPAC nucleotide codes they want. 

`make_data.sh` will create the dataset by calling `datagen.py`; generates sequences with 0 A and 0 B motifs, only A, only B, and both A and B motifs. Currently, if a motif is specified, generates up to 20 motifs. 

`datagen.py` Creates sequences with the specified number of motifs. Present parameters include:
- Nseqs = 1e4
- min_seqlen = 1e3
- max_seqlen = 1e4
- motifA = "CGACCGAACTCC"
- motifB = "ACATGCTTAGTA"

`check_data` Ensures the number of motifs desired is in the dataset matches labels (unneeded but good to test)

## Binary Classification Model [models/binary_model.py]
Identify if a specified motif is present or omitted from the data.

## Multi-Classification Model [models/multi_model.py]
Identify if only A, only B, both A/B, or no A/B motifs are present in the data.


### Miscellanous
[Random Notes]  
One possible dataset provided for training is from a paper. To avoid querying online, I downloaded it from here:
Training Data:https://raw.githubusercontent.com/abidlabs/deep-learning-genomics-primer/master/sequences.txt
Labeled Data: https://raw.githubusercontent.com/abidlabs/deep-learning-genomics-primer/master/labels.txt


The supplemental paper to accompany it is: https://www.nature.com/articles/s41588-018-0295-5
Citation: Zou, J., Huss, M., Abid, A. et al. A primer on deep learning in genomics. Nat Genet 51, 12–18 (2019). https://doi.org/10.1038/s41588-018-0295-5


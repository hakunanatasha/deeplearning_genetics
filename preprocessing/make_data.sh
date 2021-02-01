#!/usr/bin/env bash

# Number of sequences
Nseqs=10000
# Maximum number of motifs present in the data for A/B (can be 1->maxA/maxB)
# Guaranteed at LEAST 1
maxA=20
maxB=20

min_len=1000
max_len=5000

savename_noAB="../data/motif_noA_noB"
savename_A="../data/motif_A"
savename_B="../data/motif_B"
savename_AB="../data/motif_AB"

python datagen.py $Nseqs 0 0 $savename_noAB $min_len $max_len
python datagen.py $Nseqs $maxA 0 $savename_A $min_len $max_len 
python datagen.py $Nseqs 0 $maxB $savename_B $min_len $max_len
python datagen.py $Nseqs $maxA $maxB $savename_AB $min_len $max_len

python check_data.py
#!/usr/bin/env bash

# Maximum number of motifs present in the data for A/B
maxA=20
maxB=20

python datagen.py 0 0 motif_noA_noB
python datagen.py $maxA 0 motif_A
python datagen.py 0 $maxB motif_B
python datagen.py $maxA $maxB motif_AB

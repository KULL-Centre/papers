#!/bin/sh
#

python=/lindorffgrp-isilon/thomasen/software/miniconda3/bin/python3.7

for protein in TIA1 hnRNPA1 hSUMO_hnRNPA1
do

cd $protein
cp ../PW_rescaling_martini3.py .

for l in 1.00 1.02 1.04 1.06 1.08 1.10 1.12 1.14
do
$python PW_rescaling_martini3.py -i dihedrals_rubberbands_all_PRO.top -o all_PRO_lambda${l}.top -l $l

done

cd ..

done

#!/bin/sh
#

python=/storage1/thomasen/software/miniconda3/bin/python3.7

for protein in aSyn FUS ubq
do

cd $protein/two_${protein}_init
cp ../../PW_rescaling_martini3.py .

for l in 1.00 1.06 1.08
do

$python PW_rescaling_martini3.py -i all_PRO.top -o all_PRO_lambda${l}.top -l $l

done

cd ..

done


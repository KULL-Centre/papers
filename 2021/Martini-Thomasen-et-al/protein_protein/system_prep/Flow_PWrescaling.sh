#!/bin/sh
#

python=/lindorffgrp-isilon/thomasen/software/miniconda3/bin/python3.7

for protein in p15PAF htau40 FUS aSyn
do

cd two_${protein}_init
cp ../PW_rescaling_martini3.py .

for j in 1.00 1.10 1.12
do
$python PW_rescaling_martini3.py -i all_PRO.top -o all_PRO_lambda${j}.top -l $j -n 2

done

cd ..

done

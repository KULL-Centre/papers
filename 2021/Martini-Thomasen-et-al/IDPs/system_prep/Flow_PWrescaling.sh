#!/bin/sh
#

python=/storage1/thomasen/software/miniconda3/bin/python3.7

for protein in FUS A2 K25 PNt A1 CoRNID ColNT FhuA Hst52 K19 PNt Sic1 aSyn Hst5 ACTR 
do

cd $protein
cp ../PW_rescaling_martini3.py .

for l in 1.00 1.04 1.06 1.08 1.10 1.12 1.14
do

$python PW_rescaling_martini3.py -i all_PRO.top -o all_PRO_lambda${l}.top -l $l

done

cd ..

done


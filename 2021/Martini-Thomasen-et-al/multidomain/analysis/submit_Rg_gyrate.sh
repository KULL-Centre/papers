#!/bin/bash

for protein in hSUMO_hnRNPA1 hnRNPA1 TIA1
do

cd $protein

for j in 1.00 1.04 1.06 1.08 1.10 1.12 1.14
do

cd lambda_${j}
cp ../../make_Rg_gyrate.sh .
qsub make_Rg_gyrate.sh
rm \#*

cd ..

done

cd ..

done


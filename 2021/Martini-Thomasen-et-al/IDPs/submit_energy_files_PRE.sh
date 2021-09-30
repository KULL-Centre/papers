#!/bin/bash

for i in FUS A2 aSyn_PRE
do

cd $i

for j in 1.00 1.06 1.08
do

cd lambda_${j}
cp ../../make_energy_files.sh .
qsub make_energy_files.sh
cd ..

done

cd ..

done


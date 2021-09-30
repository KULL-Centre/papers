#!/bin/bash

for protein in ubq FUS aSyn
do

cd $protein

for k in $(seq 1 10)
do

cd two_${protein}_${k}

for j in 1.00 1.06 1.08
do

cd lambda_${j}
cp ../../../process_traj_pbc.sh .
qsub process_traj_pbc.sh
rm \#*
cd ..

done

cd ..

done

cd ..

done

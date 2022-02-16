#!/bin/bash

for i in FUS htau40 p15PAF aSyn ubq villin_h36
do

cd $i

for k in $(seq 1 10)
do

cd two_${i}_${k}

for j in 1.00 1.10 1.12
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

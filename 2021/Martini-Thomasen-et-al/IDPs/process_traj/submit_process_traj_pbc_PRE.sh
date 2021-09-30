#!/bin/bash

for i in A2 FUS aSyn_PRE
do

cd $i

for j in 1.00 1.06 1.08
do

cd lambda_${j}
cp ../../process_traj_pbc.sh .
qsub process_traj_pbc.sh
rm \#*
cd ..

done

cd ..

done


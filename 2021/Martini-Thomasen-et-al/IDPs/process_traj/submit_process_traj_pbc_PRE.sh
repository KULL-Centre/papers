#!/bin/bash

for i in htau40 OPN aSyn_PRE FUS A2
do

cd $i

for j in 1.00 1.10 1.12
do

cd lambda_${j}
cp ../../process_traj_pbc.sh .
qsub process_traj_pbc.sh
rm \#*
cd ..

done

cd ..

done


#!/bin/bash

t=298


cd two_villin_h36_init

for j in 1.00
do

cd lambda_${j}
cp ../../prodrun_grompp_initial.sh .
qsub prodrun_grompp_initial.sh -v temp=$t
cd ..

done

cd ..



#!/bin/bash

t=298


cd two_villin_h36_init

for j in 1.00
do

mkdir lambda_$j
cd lambda_$j
cp ../../relax_grompp.sh .
cp ../all_PRO_lambda${j}.top .
mv all_PRO_lambda${j}.top all_PRO_lambda.top
qsub relax_grompp.sh -v temp=$t
cd ..

done

cd ..



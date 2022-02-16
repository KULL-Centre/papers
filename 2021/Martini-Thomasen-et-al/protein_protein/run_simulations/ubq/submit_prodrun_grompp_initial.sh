#!/bin/bash

t=303


cd two_ubq_init

for j in 1.00
do

cd lambda_${j}
cp ../../prodrun_grompp_initial.sh .
qsub prodrun_grompp_initial.sh -v temp=$t
cd ..

done

cd ..



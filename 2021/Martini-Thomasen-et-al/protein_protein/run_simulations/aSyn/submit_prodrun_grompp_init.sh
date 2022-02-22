#!/bin/bash

t=283


cd two_aSyn_init

for j in 1.00
do

cd lambda_${j}
cp ../../prodrun_grompp_init.sh .
qsub prodrun_grompp_init.sh -v temp=$t
cd ..

done

cd ..



#!/bin/bash

t=298

cd two_FUS_init

for j in 1.00
do

cd lambda_${j}
cp ../../prodrun_grompp_init.sh .
qsub prodrun_grompp_init.sh -v temp=$t
cd ..

done

cd ..



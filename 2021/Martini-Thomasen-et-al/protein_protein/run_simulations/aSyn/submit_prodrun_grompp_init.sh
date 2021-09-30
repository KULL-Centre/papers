#!/bin/bash

t=283


cd two_aSyn_init

cd lambda_1.00
cp ../../prodrun_grompp_init.sh .
qsub prodrun_grompp_init.sh -v temp=$t
cd ..

cd ..



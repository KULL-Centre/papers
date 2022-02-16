#!/bin/bash

t=283


cd two_aSyn_init


mkdir lambda_1.00
cd lambda_1.00
cp ../../relax_grompp.sh .
cp ../all_PRO_lambda1.00.top .
mv all_PRO_lambda1.00.top all_PRO_lambda.top
qsub relax_grompp.sh -v temp=$t
cd ..


cd ..



#!/bin/bash

for i in hSUMO_hnRNPA1
do

cd $i

t=0

if [[ "$i" == "TIA1" || "$i" == "hnRNPA1" || "$i" == "hSUMO_hnRNPA1" ]]
then
        t=300
fi

echo "$i is at ${t}K"

for j in 1.00 1.02 1.04 1.06 1.08 1.10 1.12 1.14
do

mkdir lambda_$j
cd lambda_$j
cp ../../relax_grompp_hSUMO_hnRNPA1.sh .
cp ../all_PRO_lambda${j}.top .
mv all_PRO_lambda${j}.top all_PRO_lambda.top
qsub relax_grompp_hSUMO_hnRNPA1.sh -v temp=$t
cd ..

done

cd ..

done


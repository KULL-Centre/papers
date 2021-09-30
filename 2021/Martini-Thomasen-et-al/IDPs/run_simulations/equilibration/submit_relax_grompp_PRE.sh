#!/bin/bash

for i in aSyn_PRE FUS A2
do

cd ${i}

if [[ "$i" == "aSyn_PRE" ]]
then
        t=283
fi

if [[ "$i" == "FUS" || "$i" == "A2" ]]
then
        t=298
fi

echo "$i is at ${t}K"

for j in 1.00 1.06 1.08
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

done


#!/bin/bash

for i in FUS htau40 OPN aSyn_PRE A2
do

cd ${i}

t=0

if [[ "$i" == "htau40" ]]
then
        t=278
fi

if [[ "$i" == "aSyn_PRE" ]]
then
        t=283
fi

if [[ "$i" == "OPN" || "$i" == "FUS" || "$i" == "A2" ]]
then
        t=298
fi

echo "$i is at ${t}K"

for j in 1.00 1.10 1.12
do

cd lambda_${j}
cp ../../prodrun_grompp.sh .
qsub prodrun_grompp.sh -v temp=$t
cd ..

done

cd ..

done


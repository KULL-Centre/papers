#!/bin/bash

for i in aSyn_PRE FUS A2
do

cd ${i}


if [[ "$i" == "aSyn" ]]
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

cd lambda_${j}
cp ../../prodrun_grompp.sh .
qsub prodrun_grompp.sh -v temp=$t
cd ..

done

cd ..

done


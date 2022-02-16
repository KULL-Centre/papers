#!/bin/bash

for i in Hst5_PME aSyn_PME
do

t=0

if [[ "$i" == "Hst5_PME" || "$i" == "aSyn_PME" ]]
then
        t=293
fi

cd $i
echo "$i is at ${t}K"


for j in 1.00 1.04 1.06 1.08 1.10 1.12 1.14
do

cd lambda_${j}
cp ../../prodrun_grompp_PME.sh .
qsub prodrun_grompp_PME.sh -v temp=$t
cd ..

done

cd ..

done


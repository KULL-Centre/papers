#!/bin/bash

for i in TIA1 hnRNPA1 hSUMO_hnRNPA1
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

cd lambda_${j}
cp ../../prodrun_grompp.sh .
qsub prodrun_grompp.sh -v temp=$t
cd ..

done

cd ..

done


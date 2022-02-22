#!/bin/bash

for protein in htau40 OPN aSyn_PRE FUS A2
do

for l in 1.00 1.10 1.12
do 

mkdir $protein/lambda_$l/Backmapping
cd $protein/lambda_$l/Backmapping
cp ../../../Backmap.py .

echo "#!/bin/sh" > temp
echo "#PBS -W group_list=ku_10001 -A ku_10001" >> temp
echo "#PBS -N ${protein}_${l}_backmap" >> temp
cat ../../../Backmap.sh >> temp
mv temp Backmap.sh

qsub Backmap.sh

cd ../../..

done

done


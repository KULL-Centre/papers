#!/bin/bash

for i in $(seq 1 10)
do

for l in 1.00 1.10 1.12
do

#mkdir two_FUS_$i/lambda_$l/Backmapping
cd two_FUS_$i/lambda_$l/Backmapping
cp ../../../Backmap.py .

echo "#!/bin/sh" > temp
echo "#PBS -W group_list=ku_10001 -A ku_10001" >> temp
echo "#PBS -N twoFUS_${i}_${l}_backmap" >> temp
cat ../../../Backmap.sh >> temp
mv temp Backmap.sh

qsub Backmap.sh

cd ../../..

done

done


#!/bin/bash

for i in $(seq 1 10)
do

cd two_ubq_${i}

for j in 1.00 1.10 1.12
do

cd lambda_$j

cp ../../energy.sh .
qsub energy.sh
rm \#*

cd ..

done

cd ..

done
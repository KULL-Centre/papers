#!/bin/bash

for i in FUS A2 aSyn_PRE
do

cd $i

for j in 1.00 1.06 1.08
do

cd lambda_${j}
cp ../../make_PI_mindist.sh .
qsub make_PI_mindist.sh
rm \#*

cd ..

done

cd ..

done


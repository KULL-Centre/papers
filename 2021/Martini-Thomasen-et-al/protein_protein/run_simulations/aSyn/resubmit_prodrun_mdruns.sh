#!/bin/bash 
 
for i in $(seq 1 10)
do 
 
cd two_aSyn_$i 
 
for j in 1.00 1.06 1.08
do 
 
cd lambda_${j} 
 
qsub prodrun_mdrun.sh
cd .. 
 
done 
 
cd .. 
 
done 

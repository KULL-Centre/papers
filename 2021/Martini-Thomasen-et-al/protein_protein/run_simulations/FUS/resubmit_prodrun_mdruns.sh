#!/bin/bash 
 
for i in $(seq 1 10)
do 
 
cd two_FUS_$i 
 
for j in 1.00 1.10 1.12 
do 
 
cd lambda_${j} 
qsub prodrun_mdrun.sh
cd .. 
 
done 
 
cd .. 
 
done 

#!/bin/bash 
 
for i in aSyn_PRE FUS A2
do 
 
cd $i 
 
for j in 1.00 1.06 1.08
do 
 
cd lambda_${j} 
 
qsub prodrun_mdrun.sh
cd .. 
 
done 
 
cd .. 
 
done 

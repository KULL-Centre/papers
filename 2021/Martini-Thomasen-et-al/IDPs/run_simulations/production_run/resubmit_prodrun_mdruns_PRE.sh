#!/bin/bash 
 
for i in htau40 OPN aSyn_PRE FUS A2
do 
 
cd $i 
 
for j in 1.00 1.10 1.12
do 
 
cd lambda_${j} 
 
qsub prodrun_mdrun.sh
cd .. 
 
done 
 
cd .. 
 
done 

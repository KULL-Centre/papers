#!/bin/bash 
 
for i in hSUMO_hnRNPA1 hnRNPA1 TIA1
do 
 
cd $i 
 
for j in 1.00 1.04 1.06 1.08 1.10 1.12 1.14
do 
 
cd lambda_${j} 
 
qsub prodrun_mdrun.sh
cd .. 
 
done 
 
cd .. 
 
done 

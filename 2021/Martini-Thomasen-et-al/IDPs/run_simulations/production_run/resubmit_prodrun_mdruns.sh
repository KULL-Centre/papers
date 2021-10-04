#!/bin/bash 
 
for i in K25 A1 ColNT FhuA Hst52 K19 Sic1 aSyn Hst5 ACTR CoRNID PNt
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

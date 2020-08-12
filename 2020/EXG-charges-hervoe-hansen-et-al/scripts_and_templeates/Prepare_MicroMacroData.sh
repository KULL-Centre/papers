#!/bin/bash
#
# Collecting conditional probabilities of the various protonation state combinations

echo "% pH      PD         PP         DP         DD         Total       Site1       Site2" >data_Ullmann.dat
for pH in 2.00 2.50 3.00 3.50 4.00 4.50 5.00 5.50 6.00 6.50 7.00 7.50 8.00 8.50 9.00 9.50
do

awk 'NR>1 && NR<=5' pH_${pH}_conditional_prob.dat |  tr '\n' ' | ' | sed -e 's/|$/\n/' >a
printf $pH >>data_Ullmann.dat
#2.00               39:P,66:D   0.005161               39:P,66:P   0.990582               39:D,66:P   0.004250               39:D,66:D   0.000007 
cat a | awk '{printf(" %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f \n",$2,$4,$6,$8,1*$2+2*$4+1*$6+0*$8,$2+$4,$4+$6)}' >> data_Ullmann.dat

done

cat data_Ullmann.dat

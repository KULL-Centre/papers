#!/bin/sh
# Stefan H. Hansen & Yong Wang

cpin=../EXG.cpin

rm reordered_cpouts.pH_* pH*
echo cphstats --fix-remd reordered_cpouts ../EXG_explicit_pHREMD/pH*/EXG.pHREMD.cpout 
cphstats --fix-remd reordered_cpouts ../EXG_explicit_pHREMD/pH*/EXG.pHREMD.cpout

echo "# 1    2" >titration_curves.dat
echo "# pH   Residue" >>titration_curves.dat
for pH in 0.50 1.00 1.50 2.00 2.50 3.00 3.50 4.00 4.50 5.00 5.50 6.00 6.50 7.00 7.50 8.00
do
i=`expr $i + 1`
dir="MD_pH${pH}"
   	cphstats -i $cpin reordered_cpouts.pH_${pH} -o pH${pH}_calcpka.dat --population pH${pH}_populations.dat

	cphstats -i $cpin reordered_cpouts.pH_${pH} -n 10000 --cumulative --cumulative-out pH_${pH}_cumulative.dat
	cphstats -i $cpin reordered_cpouts.pH_${pH} -n 10000 -r 100000 -R pH_${pH}_runningavg.dat 

	echo $pH | awk '{printf("%4.2f ",$1)}'>> titration_curves.dat
	cat pH${pH}_calcpka.dat | awk '/Frac Prot/ {printf("%5.3f ",1-$10)}' >>titration_curves.dat
	echo $pH | awk '{printf(" \n",$1)}'>> titration_curves.dat
done


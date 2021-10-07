#!/bin/bash
#PBS -W group_list=ku_10001 -A ku_10001
#PBS -l nodes=1:ppn=1
#PBS -l walltime=200:00:00

#cd $PBS_O_WORKDIR

p=$1

dmax=`more DMAX_$p | awk 'BEGIN{c=0}{c=c+$1}END{print c/NR}'`
cp ../SAXS/$p.dat Data.dat
echo -e "Data.dat\n\n\n$dmax\n\n\n\n\n\n\n\n\n\n" > inputfile.d
./bift < inputfile.d

mv rescale.d $p.dat

for l in `more ../LAMBDAS`
do
	python /home/projects/ku_10001/people/frapes/Martini3_PW/backmapper/saxs_chi.py ../$p/lambda_$l/calc_saxs.dat $p.dat bift_${p}_l$l.dat
done

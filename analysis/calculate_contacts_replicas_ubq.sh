#!/bin/sh
#PBS -W group_list=ku_10001 -A ku_10001
#PBS -N ubq_contacts
#PBS -l nodes=1:ppn=40:thinnode
#PBS -l walltime=48:00:00
#PBS -l mem=130gb
# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes
# Load all required modules for the job
module load tools
module load cuda/toolkit/10.2.89 openmpi/gcc/64/1.10.2 gcc/9.3.0
gmx=/home/projects/ku_10001/apps/GMX20203/bin/gmx_mpi

for j in $(seq 1 10)
do

cd two_ubq_$j

for i in 1.00 1.10 1.12
do

cd lambda_${i}

mkdir data

$gmx make_ndx -f prodrun.tpr -o data/ubq1_ubq2_lambda${i}.ndx <<EOF
a 1-166
a 167-332
q
EOF
 
$gmx mindist -f prodrun_nopbc.xtc -s prodrun.tpr -n data/ubq1_ubq2_lambda${i}.ndx -od data/ubq1_ubq2_mindist_lambda${i}.xvg -on data/ubq1_ubq2_numcont_lambda${i}.xvg -tu us -d 0.8 <<EOF 
a_1-166
a_167-332
EOF

cd ..

done 

cd ..

done

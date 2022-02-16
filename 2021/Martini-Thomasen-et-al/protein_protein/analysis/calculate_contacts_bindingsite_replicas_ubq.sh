#!/bin/sh
#PBS -W group_list=ku_10001 -A ku_10001
#PBS -N ubq_contacts_bindingsite
#PBS -l nodes=1:ppn=40
#PBS -l walltime=48:00:00
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


#Binding site residue 8, 13, 44, 45, 46, 49, 67, 68, 70, 71, and 73

for j in $(seq 1 10)
do

cd two_ubq_$j


for i in 1.00 1.10 1.12
do

cd lambda_$i

$gmx make_ndx -f prodrun.tpr -o data/ubq1_ubq2_lambda${i}_bindingsite.ndx <<EOF
ri 8 | ri 13 | ri 44-46 | ri 49 | ri 67-68 | ri 70-71 | ri 73
ri 84 | ri 89 | ri 120-122 | ri 125 | ri 143-144 | ri 146-147 | ri 149
q

EOF
 
$gmx mindist -f prodrun_nopbc.xtc -s prodrun.tpr -n data/ubq1_ubq2_lambda${i}_bindingsite.ndx -od data/ubq1_ubq2_mindist_lambda${i}_bindingsite.xvg -on data/ubq1_ubq2_numcont_lambda${i}_bindingsite.xvg -tu us -d 0.8 <<EOF 
r_8_r_13_r_44-46_r_49_r_67-68_r_70-71_r_73
r_84_r_89_r_120-122_r_125_r_143-144_r_146-147_r_149
EOF

cd ..

done 

cd ..

done

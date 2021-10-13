#!/bin/bash
#PBS -W group_list=ku_10001 -A ku_10001
#PBS -l nodes=1:ppn=10
#PBS -l walltime=200:00:00

cd $PBS_O_WORKDIR

module purge
module load tools
module load cuda/toolkit/10.2.89 openmpi/gcc/64/1.10.2
module load gcc/9.3.0
source /home/projects/ku_10001/apps/gromacs-2018.1_AVX2_256_12JUN2018/bin/GMXRC
scripts=/home/projects/ku_10001/people/frapes/Martini3_PW/backmapper

protein=$1
lambda=$2

#rm -r $protein

mkdir $protein
mkdir $protein/lambda_$lambda
cd $protein/lambda_$lambda

#EXTRACT FRAMES
mkdir CG_ENSEMBLE
python $scripts/extract.py /home/projects/ku_10001/people/fretho/MARTINI_PWrescaling/IDP/$protein/lambda_$lambda
/prodrun_nopbc.xtc /home/projects/ku_10001/people/fretho/MARTINI_PWrescaling/IDP/$protein/PRO_CG.pdb 15000

#BACKMAPPING AND SAXS
gmx_mpi pdb2gmx -f /home/projects/ku_10001/people/fretho/MARTINI_PWrescaling/IDP/$protein/min_AA.pdb -o AA.gro -
ignh -water none -ff charmm27 -quiet
mkdir AA_ENSEMBLE
ens_size=`find CG_ENSEMBLE/ -name "CG_frame*pdb" | wc -l`
exp_saxs=/home/projects/ku_10001/people/frapes/Martini3_PW/SAXS/$protein.dat
i=0
while [ $i -lt $ens_size ]
do
	bash $scripts/initram-v5_200.sh -f CG_ENSEMBLE/CG_frame$i.pdb -p topol.top -to charmm36 >/dev/null
	gmx_mpi editconf -f backmapped.gro -o AA_ENSEMBLE/AA_frame$i.pdb -quiet
	/home/projects/ku_10001/apps/Pepsi-SAXS AA_ENSEMBLE/AA_frame$i.pdb $exp_saxs -o saxs.dat -cst --cstFacto
r 0 --I0 1.0 --dro 1.0 --r0_min_factor 1.025 --r0_max_factor 1.025 --r0_N 1
	echo `cat saxs.dat | awk '{if (NF==4) printf $4" "}END{printf "\n"}'` >> calc_saxs.dat
	i=$(( $i + 1 ))
done

python $scripts/saxs_chi.py /home/projects/ku_10001/people/frapes/Martini3_PW/$protein/lambda_$lambda/calc_saxs.dat /home/projects/ku_10001/people/frapes/Martini3_PW/SAXS/$protein.dat

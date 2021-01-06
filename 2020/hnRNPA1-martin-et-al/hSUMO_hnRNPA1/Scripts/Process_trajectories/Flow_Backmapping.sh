#!/bin/bash
#SBATCH --job-name=hnRNPA1_backmapping
#SBATCH --time=144:00:00
#SBATCH --partition=sbinlab
#SBATCH --nodes=1 --ntasks-per-node=64 --cpus-per-task=1
#SBATCH --mem=8Gb

export PATH=/groups/sbinlab/software/openmpi-2.1.0/bin:$PATH
export LD_LIBRARY_PATH=/groups/sbinlab/software/openmpi-2.1.0/lib:$LD_LIBRARY_PATH
export PATH=/lustre/hpc/sbinlab/software/miniconda3/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/hpc/sbinlab/software/miniconda3/lib:/groups/sbinlab/wyong/usr/local/PLUMED253/lib
export PLUMED_KERNEL=/groups/sbinlab/wyong/usr/local/PLUMED253/lib/libplumedKernel.so


gmx=/groups/sbinlab/wyong/usr/local/GMX2019.4/bin/gmx_mpi
BM=/groups/sbinlab/thomasen/scripts/backward-v5/initram-v5_200.sh
AApdb=/groups/sbinlab/thomasen/hnRNPA1_ionicstrength_productionrun_theta107/hnRNPA1_sumo_homology_model_1.pdb
FF=charmm27


for t in $1
do
	cd ${t}mM
	mkdir Backmapping_frames
	cd Backmapping_frames
	echo 1 | $gmx trjconv -f ../hnRNPA1_${t}mM_theta1.07.xtc -o frame.pdb -pbc whole -sep -s ../hnRNPA1_${t}mM_theta1.07.tpr

	Nframe=`ls -1q frame*.pdb | wc -l`

	$gmx pdb2gmx -f $AApdb -o AA.gro -ignh -water none -ff $FF -quiet

	i=0
	while [ $i -lt $Nframe ] 
	do
		name=frame${i}
		CGgro=${name}.pdb
		echo $name
		bash $BM -f $CGgro -p topol.top -to charmm36 >/dev/null
		$gmx editconf -f backmapped.gro -o AA_${name}.pdb -quiet
		i=$(( $i + 1 ))
	done


	rm #*
	rm backmapped.*


	cd ../..
done

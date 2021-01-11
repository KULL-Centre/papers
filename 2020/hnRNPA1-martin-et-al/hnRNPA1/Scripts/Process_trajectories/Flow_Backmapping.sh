#!/bin/bash


gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi
BM=/lindorffgrp-isilon/thomasen/scripts/backward-v5/initram-v5_200.sh
AApdb=/lindorffgrp-isilon/thomasen/hnRNPA1/hnRNPA1_ionicstrength_notag/hnRNPA1_noTag_homology_model_1_folded_LCD.pdb
FF=charmm27


for t in $1 
do
	cd ${t}mM
	mkdir Backmapping_frames
	cd Backmapping_frames
	echo 1 | $gmx trjconv -f ../hnRNPA1_notag_${t}mM.xtc -o frame.pdb -pbc whole -sep -s ../hnRNPA1_notag_${t}mM.tpr

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

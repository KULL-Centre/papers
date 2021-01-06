#!/bin/bash
#

gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi  


for t in 50 100 200 300 400 500
do

cd ${t}mM
	
$gmx trjconv -f hnRNPA1_${t}mM_theta1.07.xtc -s hnRNPA1_${t}mM_theta1.07.tpr -pbc mol -o tmp.xtc <<EOF
	1
	1
EOF

$gmx trjconv -f tmp.xtc -s hnRNPA1_${t}mM_theta1.07.tpr -fit rot+trans -o ./../../ionicstrength_2_trajectories_processed/hnRNPA1_${t}mM_theta1.07_processed.xtc <<EOF
	1
	1
EOF

rm tmp.xtc

cd ..

done

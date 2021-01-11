#!/bin/bash
#

gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi 

for t in 50 150 250 400 1000
do

cd ${t}mM
	
$gmx trjconv -f hnRNPA1_notag_${t}mM.xtc -s hnRNPA1_notag_${t}mM*.tpr -pbc mol -o tmp.xtc <<EOF
	1
	1
EOF

$gmx trjconv -f tmp.xtc -s hnRNPA1_notag_${t}mM*.tpr -fit rot+trans -o ../Processed_traj/hnRNPA1_${t}mM_processed.xtc <<EOF
	1
	1
EOF

rm tmp.xtc


cd ..

done

#!/bin/bash

gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi

for i in 1.00 1.04 1.06 1.07 1.08 1.09 1.10 1.12
do

	inputfile=hnRNPA1*theta${i}*

	#protein Rg
	$gmx gyrate -f Processed_traj/${inputfile}.xtc -s theta${i}/${inputfile}.tpr -o data/Rg_gyrate_theta${i}.xvg -b 1000000 <<EOF
	1
EOF

done

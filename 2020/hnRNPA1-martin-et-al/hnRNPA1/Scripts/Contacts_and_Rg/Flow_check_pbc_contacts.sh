#!/bin/bash
#

gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi

for i in 50 150 250 400 1000
do
	$gmx mindist -f hnRNPA1_${i}mM_processed.xtc -s ../${i}mM/hnRNPA1_notag_${i}mM_gpu.tpr -od check_pbc_contacts_${i}mM.xvg -tu us -pi <<EOF
	1
EOF

done

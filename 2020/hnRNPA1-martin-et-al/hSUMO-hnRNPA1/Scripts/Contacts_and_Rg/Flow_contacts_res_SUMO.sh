#!/bin/bash
#
#SUMO: resi MET1 to ILE117 

gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi

mkdir ../data/resindxfiles_SUMO
mkdir ../data/resicontacts_SUMO
mkdir ../data/residistance_SUMO

contactscutoff=0.5
timeunit=us

for t in 50 100 200 300 400 500 
do

#define file names
inputfile=hnRNPA1*_${t}mM*

mkdir ../data/resicontacts_SUMO/${t}mM
mkdir ../data/residistance_SUMO/${t}mM

for i in {1..117}

do

$gmx make_ndx -f $inputfile.tpr -o ../data/resindxfiles_SUMO/resindex${i}_${t}mM.ndx <<EOF
r $i
r 314-433
q
EOF

$gmx mindist -f ${inputfile}_processed.xtc -s ../hnRNPA1_ionicstrength_2_productionrun/${t}mM/$inputfile.tpr -n ../data/resindxfiles_SUMO/resindex${i}_${t}mM.ndx -od ../data/residistance_SUMO/${t}mM/LCD_HisSUMO_${t}mM_resi${i}.xvg -on ../data/resicontacts_SUMO/${t}mM/LCD_HisSUMO_${t}mM_resi${i}.xvg -d $contactscutoff -tu $timeunit -b 1 -group <<EOF
r_$i
r_314-433
EOF

done

done

rm -r ../data/resindxfiles_SUMO


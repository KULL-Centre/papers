#!/bin/bash
#


gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi

contactscutoff=0.5
timeunit=us

mkdir ../data/resindxfiles_LCDres
mkdir ../data/resicontacts_LCDres
mkdir ../data/residistance_LCDres

for t in 50 100 200 300 400 500
do

#define file names
inputfile=hnRNPA1*_${t}mM*   ;name of .xtc and .tpr files

mkdir ../data/resicontacts_LCDres/${t}mM
mkdir ../data/residistance_LCDres/${t}mM

for i in {314..433}

do

$gmx make_ndx -f ../hnRNPA1_ionicstrength_2_productionrun/${t}mM/$inputfile.tpr -o ../data/resindxfiles_LCDres/resindex${i}_${t}mM.ndx <<EOF
r $i
r 133-216 | r 224-303
q
EOF

$gmx mindist -f ${inputfile}_processed.xtc -s ../hnRNPA1_ionicstrength_2_productionrun/${t}mM/$inputfile.tpr -n ../data/resindxfiles_LCDres/resindex${i}_${t}mM.ndx -od ../data/residistance_LCDres/${t}mM/LCD_RRM_${t}mM_resi${i}.xvg -on ../data/resicontacts_LCDres/${t}mM/LCD_RRM_${t}mM_resi${i}.xvg -d $contactscutoff -tu $timeunit -b 1 -group <<EOF
r_$i
r_133-216_r_224-303
EOF

done

done

rm -r ../data/resindxfiles_LCDres

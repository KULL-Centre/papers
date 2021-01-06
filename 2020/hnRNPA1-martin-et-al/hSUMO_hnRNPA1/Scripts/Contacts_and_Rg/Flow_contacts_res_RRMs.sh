#!/bin/bash
#

gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi

contactscutoff=0.5
timeunit=us

mkdir ../data/resindxfiles_LCD_RRMsres
mkdir ../data/resicontacts_LCD_RRMsres
mkdir ../data/residistance_LCD_RRMsres

for t in 50 100 200 300 400 500 
do

#define file names
inputfile=hnRNPA1*_${t}mM*    ;name of .xtc and .tpr files

mkdir ../data/resicontacts_LCD_RRMsres/${t}mM
mkdir ../data/residistance_LCD_RRMsres/${t}mM

for i in {133..303}

do

$gmx make_ndx -f ../hnRNPA1_ionicstrength_2_productionrun/${t}mM/$inputfile.tpr -o ../data/resindxfiles_LCD_RRMsres/resindex${i}_${t}mM.ndx <<EOF
r $i
r 314-433
q
EOF

$gmx mindist -f ${inputfile}_processed.xtc -s ../hnRNPA1_ionicstrength_2_productionrun/${t}mM/$inputfile.tpr -n ../data/resindxfiles_LCD_RRMsres/resindex${i}_${t}mM.ndx -od ../data/residistance_LCD_RRMsres/${t}mM/LCD_RRMs_${t}mM_resi${i}.xvg -on ../data/resicontacts_LCD_RRMsres/${t}mM/LCD_RRMs_${t}mM_resi${i}.xvg -d $contactscutoff -tu $timeunit -b 1 -group <<EOF
r_$i
r_314-433
EOF

done

done

rm -r ../data/resindxfiles_LCD_RRMsres

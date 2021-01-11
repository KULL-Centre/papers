#!/bin/bash
#


gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi



contactscutoff=0.5
timeunit=us

mv ../data/resicontacts_LCDRRM ../data/resicontacts_LCDRRM_nogroup
mv ../data/residistance_LCDRRM ../data/residistance_LCDRRM_nogroup

mkdir ../data/resindxfiles_LCDRRM
mkdir ../data/resicontacts_LCDRRM
mkdir ../data/residistance_LCDRRM

for t in 50 150 250 400 1000
do

#define file names
inputfile=hnRNPA1*_${t}mM*   ;name of .xtc and .tpr files

mkdir ../data/resicontacts_LCDRRM/${t}mM
mkdir ../data/residistance_LCDRRM/${t}mM

for i in {133..303}

do

$gmx make_ndx -f ../${t}mM/$inputfile.tpr -o ../data/resindxfiles_LCDRRM/resindex${i}_${t}mM.ndx <<EOF
r $i
r 314-433
q
EOF

$gmx mindist -f ${inputfile}_processed.xtc -s ../${t}mM/$inputfile.tpr -n ../data/resindxfiles_LCDRRM/resindex${i}_${t}mM.ndx -od ../data/residistance_LCDRRM/${t}mM/LCD_RRM_${t}mM_resi${i}.xvg -on ../data/resicontacts_LCDRRM/${t}mM/LCD_RRM_${t}mM_resi${i}.xvg -d $contactscutoff -tu $timeunit -b 1 -group <<EOF
r_$i
r_314-433
EOF

done

done

rm -r ../data/resindxfiles_LCDRRM


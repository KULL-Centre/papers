#!/bin/bash
#

gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi


contactscutoff=0.5

timeunit=us

for i in 50 150 250 400 1000
do

inputfile=hnRNPA1*_${i}mM*

#protein Rg
$gmx gyrate -f ${inputfile}.xtc -s ../${i}mM/${inputfile}.tpr -o ../data/Rg_gyrate_${i}mM.xvg -b 1000000 <<EOF
1
EOF

#make groups for calculating contacts 
$gmx make_ndx -f ../${i}mM/$inputfile.tpr -o ../${i}mM/LCDRRMdomains.ndx <<EOF
r 133-303
r 314-433
q
EOF

#minimum distance and number of contacts between defined groups
$gmx mindist -f ${inputfile}_processed.xtc -s ../${i}mM/$inputfile.gro -n ../${i}mM/LCDRRMdomains.ndx -od ../data/LCDRRM_mindist_minusstart_${i}mM.xvg -on ../data/LCDRRM_contacts_minusstart_${i}mM.xvg -d $contactscutoff -tu $timeunit -b 1 -group <<EOF
r_133-303
r_314-433
EOF

done


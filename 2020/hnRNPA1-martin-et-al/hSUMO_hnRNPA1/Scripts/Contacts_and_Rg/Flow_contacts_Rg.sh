#!/bin/bash
#

gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi


#define contact groups and cutoff distance
contactscutoff=0.5

#define unit of time
timeunit=us


for t in 50 100 200 300 400 500 
do

#define file names
inputfile=hnRNPA1_${t}mM_theta1.07   ;name of .xtc and .tpr files

#mkdir xvg_files_$inputfile

#protein Rg
$gmx gyrate -f ${inputfile}_processed.xtc -s $inputfile.tpr -o /lindorffgrp-isilon/thomasen/hnRNPA1/hnRNPA1_ionicstrength_2/data/Rg_minusstart_$inputfile.xvg -b 1000000 <<EOF
1
EOF 

#minimum distance and contacts with periodic mirror image
$gmx mindist -f ${inputfile}_processed.xtc -s $inputfile.tpr -pi -od /lindorffgrp-isilon/thomasen/hnRNPA1/hnRNPA1_ionicstrength_2/data/mindist_pi_minusstart_$inputfile.xvg -tu $timeunit <<EOF
1
EOF

#make groups for calculating contacts 
$gmx make_ndx -f $inputfile.tpr -o SUMOandRRMdomains.ndx <<EOF
r 1-117 | r 133-303
r 314-433
q
EOF

#minimum distance and number of contacts between defined groups
$gmx mindist -f ${inputfile}_processed.xtc -s $inputfile.tpr -n SUMOandRRMdomains.ndx -od /lindorffgrp-isilon/thomasen/hnRNPA1/hnRNPA1_ionicstrength_2/data/mindist_SUMORRM_LCD_minusstart_$inputfile.xvg -on /lindorffgrp-isilon/thomasen/hnRNPA1/hnRNPA1_ionicstrength_2/data/numcont_SUMORRM_LCD_minusstart_$inputfile.xvg -d $contactscutoff -tu $timeunit -b 1 -group <<EOF
r_1-117_r_133-303
r_314-433
EOF

#make groups for calculating contacts 
$gmx make_ndx -f $inputfile.tpr -o RRMdomains.ndx <<EOF
r 133-303
r 314-433
q
EOF

#minimum distance and number of contacts between defined groups
$gmx mindist -f ${inputfile}_processed.xtc -s $inputfile.tpr -n RRMdomains.ndx -od /lindorffgrp-isilon/thomasen/hnRNPA1/hnRNPA1_ionicstrength_2/data/mindist_allRRM_LCD_minusstart_$inputfile.xvg -on /lindorffgrp-isilon/thomasen/hnRNPA1/hnRNPA1_ionicstrength_2/data/numcont_allRRM_LCD_minusstart_$inputfile.xvg -d $contactscutoff -tu $timeunit -b 1 -group <<EOF
r_133-303
r_314-433
EOF

$gmx make_ndx -f $inputfile.tpr -o LCD_HisSUMO_index_$inputfile.ndx <<EOF
ri 314-433
ri 1-117
q
EOF

$gmx mindist -f ${inputfile}_processed.xtc -s $inputfile.tpr -n LCD_HisSUMO_index_$inputfile.ndx -od /lindorffgrp-isilon/thomasen/hnRNPA1/hnRNPA1_ionicstrength_2/data/mindist_LCD_HisSUMO_$inputfile.xvg -on /lindorffgrp-isilon/thomasen/hnRNPA1/hnRNPA1_ionicstrength_2/data/numcont_LCD_HisSUMO_$inputfile.xvg -d $contactscutoff -tu $timeunit -b 1 -group <<EOF
r_1-117
r_314-433
EOF



done

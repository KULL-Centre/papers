#!/bin/bash
#

gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi

contactscutoff=0.5

#define unit of time
timeunit=us


for t in 1.0 1.02 1.04 1.05 1.06 1.07 1.08 1.10 
do

#define file names
inputfile=hnRNPA1_theta${t}   ;name of .xtc and .tpr files
outputdir=/lindorffgrp-isilon/thomasen/hnRNPA1/hnRNPA1_PWrescaling_3/data


#mkdir xvg_files_$inputfile

#protein Rg
$gmx gyrate -f ${inputfile}_processed.xtc -s $inputfile.tpr -o ${outputdir}/Rg_minusstart_$inputfile.xvg -b 1000000 <<EOF
1
EOF 


#minimum distance and contacts with periodic mirror image
$gmx mindist -f ${inputfile}_processed.xtc -s $inputfile.tpr -pi -od ${outputdir}/mindist_pi_minusstart_$inputfile.xvg -tu $timeunit <<EOF
1
EOF

done

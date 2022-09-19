#!/bin/bash
#
#Script to coarse-grain and run energy minimization


export PATH="/lindorffgrp-isilon/thomasen/software/miniconda3/bin:$PATH"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lindorffgrp-isilon/wyong/software/openmpi401/lib

gmx=/lindorffgrp-isilon/wyong/software/GMX20194/bin/gmx_mpi
python=/lindorffgrp-isilon/thomasen/software/miniconda3/bin/python3.7

martinize=/lindorffgrp-isilon/thomasen/software/miniconda3/bin/martinize2

#wget http://cgmartini.nl/images/tools/insane/insane.py
insane=insane.py

minmdp=minimization.mdp
FF=martini3001
ffdir=/storage1/thomasen/software/force-fields/Martini/martini_v300

dssp=/lindorffgrp-isilon/thomasen/software/miniconda3/bin/mkdssp

export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

#Loop through IDPs
for i in A2 aSyn_PRE K25 A1 CoRNID ColNT FhuA Hst52 K19 PNt Sic1 aSyn Hst5 ACTR OPN htau40 FUS
do

if [[ "$i" == "A2" ]]
then
	salt=0.005
fi

if [[ "$i" == "A1" ]]
then
        salt=0.05
fi

if [[ "$i" == "htau40" ]]
then
	salt=0.1
fi

if [[ "$i" == "aSyn_PRE" ]]
then
	salt=0.125
fi

if [[ "$i" == "FhuA" || "$i" == "Hst52" || "$i" == "K19" || "$i" == "K25" || "$i" == "PNt" || "$i" == "Hst5" || "$i" == "FUS" || "$i" == "OPN" ]]
then
        salt=0.15
fi

if [[ "$i" == "CoRNID" || "$i" == "Sic1" || "$i" == "ACTR" || "$i" == "aSyn" ]]
then
        salt=0.2
fi

if [[ "$i" == "ColNT" ]]
then
        salt=0.4
fi

echo "$i has ${salt}M salt"

#Starting structure
pdb=min_AA.pdb

#Make directory, cp .pdb and insane.py file there and cd there
dir=${i}
cp $insane $dir
cd $dir
$gmx editconf -f min_AA.gro -o $pdb 

#Martinize
$python $martinize -f $pdb -o PRO_topol.top -x PRO_CG.pdb -ff $FF -ff-dir $ffdir/martini_v3.0.0_proteins/force_fields/ -map-dir $ffdir/martini_v3.0.0_proteins/mappings/ 

#Put protein in box
$gmx editconf -f PRO_CG.pdb -o PRO_CG.gro -bt dodecahedron -d 4.0 <<EOF
1
EOF

#Solvate using insane.py
python2.7 $insane -f PRO_CG.gro -o PRO_SOL_IONS.gro -pbc keep -salt $salt -sol W -center -p PRO_topol_SOL_IONS.top

#The next few blocks modify the toplogy file and molecule_0.itp file:

#Remove #include include martini.itp and substitute ion names in topology file
perl -pi -e's/#include "martini.itp"//g' PRO_topol_SOL_IONS.top
perl -pi -e's/NA\+/NA/g' PRO_topol_SOL_IONS.top
perl -pi -e's/CL-/CL/g' PRO_topol_SOL_IONS.top

#Rename molecule_0.itp to PRO.itp and rename "molecule_0" as "Protein" in PRO.itp file
mv molecule_0.itp PRO.itp
perl -pi -e's/molecule_0/Protein/g' PRO.itp

#Add "#include .itp" lines to PRO_topol_SOL_IONS.top
cat <<EOF > others.top
#include "$ffdir/martini_v3.0.0.itp"
#include "PRO.itp"
#include "$ffdir/martini_v3.0.0_ions_v1.itp"
#include "$ffdir/martini_v3.0.0_solvents_v1.itp"
EOF
cat others.top PRO_topol_SOL_IONS.top >a
mv a PRO_topol_SOL_IONS.top

#Run energy minimization
$gmx grompp -f ../$minmdp -p PRO_topol_SOL_IONS.top -c PRO_SOL_IONS.gro -o min.tpr -pp all_PRO.top -maxwarn 3 -r PRO_SOL_IONS.gro
nohup $gmx mdrun -deffnm min -v -ntomp 1 &

rm $insane

cd ..

done

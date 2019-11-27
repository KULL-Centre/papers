#!/bin/bash
#
# Flow_Backmapping.sh

BM=/lindorffgrp-isilon/andreas/prj_TIA1/TIA1_Martini3_seq2/test_backmapping/backward-v5/initram-v5_200.sh
AApdb=/lindorffgrp-isilon/andreas/TIA1_Martini3_seq2/TIA1_model2.pdb
FF=charmm27

for t in {1.00,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,1.10}
do
    mkdir theta_2_$t/Backmapping
    cd theta_2_$t/Backmapping
    echo 1 | $gmx trjconv -f ../md2.xtc -o frame.pdb -pbc whole -skip 10000 -sep -s ../md2.tpr -quiet #about 100 frames, 20 min
    
    Nframe=`ls -1q frame*.pdb | wc -l`

    $gmx pdb2gmx -f $AApdb -o AA.gro -ignh -water none -ff $FF -quiet

    i=0
    while [ $i -lt $Nframe ] 
    do
        name=frame${i}
        CGgro=${name}.pdb
        echo $name
        bash $BM -f $CGgro -p topol.top -to charmm36 >/dev/null
        $gmx editconf -f backmapped.gro -o AA_${name}.pdb -quiet
        i=$(( $i + 1 ))
    done

    rm #*
    rm backmapped.*

    cd ../..
done

#!/bin/bash

for c in 50 150 250 400 1000
do
program=run_pep_stage1_v2.py
    folder=${c}mM/stage
    rm -r $folder
    if [ ! -d "$folder" ]
    then
        if [ ! -d "${c}mM" ]
	then
            mkdir ${c}mM
	fi
	mkdir $folder
    fi
    cp $program $folder
    cd $folder
   #rm nohup.out
    nohup /storage1/thomasen/software/miniconda3/bin/python3 $program FLA1_${c}mMNaCl &
    cd ../../
done

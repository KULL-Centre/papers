#!/bin/bash

for i in 50 100 200 300 400 500
do
program=run_pep_stage1_v2.py
    folder=${i}mM/stage_reduced
    if [ ! -d "$folder" ]
    then
        if [ ! -d "${i}mM" ]
	then
            mkdir ${i}mM
	fi
	mkdir $folder
    fi
    cp $program $folder
    cd $folder
    rm nohup.out
    nohup /storage1/thomasen/software/miniconda3/bin/python3 $program $i &
    cd ../../
done

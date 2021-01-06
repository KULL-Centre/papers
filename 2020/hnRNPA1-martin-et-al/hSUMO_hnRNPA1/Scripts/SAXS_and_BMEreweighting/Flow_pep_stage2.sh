#!/bin/bash

for i in 100 200 300 400 500
do
    folder=${i}mM/stage_reduced/
    program=run_pep_stage2_v4.py
    cp $program $folder
    cd $folder
    rm nohup.out 
    nohup /storage1/thomasen/software/miniconda3/bin/python3 $program $i &
    cd ../../
done

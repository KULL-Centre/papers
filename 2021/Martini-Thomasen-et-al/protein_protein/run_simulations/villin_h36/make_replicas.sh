#!/bin/bash

python=/home/projects/ku_10001/people/fretho/miniconda3/bin/python3.8

for i in $(seq 1 10)
do

mkdir two_villin_h36_$i
cp two_villin_h36_init/all_PRO_lambda* two_villin_h36_$i

done 

$python make_replica_startstructure.py

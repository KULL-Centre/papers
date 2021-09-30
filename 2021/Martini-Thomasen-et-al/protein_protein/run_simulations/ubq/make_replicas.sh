#!/bin/bash

python=/home/projects/ku_10001/people/fretho/miniconda3/bin/python3.8

for i in $(seq 1 10)
do

mkdir two_ubq_$i
cp two_ubq_init/all_PRO_lambda* two_ubq_$i

done 

$python make_replica_startstructure.py

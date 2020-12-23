#!/bin/bash

program=run_BME_stage2.py
# j=0 for SAXS
# j=1 for SANS 0%
# j=2 for SANS 42%
# j=3 for SANS 70%
# j=4 for all
# j>4 for other combinations

#for j in {0,1,2,3,4}
for j in 0
do
    for i in 50 100 200 300 400 500 
    do
        folder=${lambda}mM/stage_reduced/.
        cp bme_reweight.py $program $folder
        cd $folder 
        /storage1/thomasen/software/miniconda2/bin/python2 $program $j 
        cd ../../
    done
done

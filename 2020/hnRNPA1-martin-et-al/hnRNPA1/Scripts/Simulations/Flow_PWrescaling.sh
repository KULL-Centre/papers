#!/bin/sh
#

export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

for i in 50 150 250 400 1000
do

theta=1.07
PWtop=all_PRO_theta${theta}_${i}mM.top
perl PWrescaling_Martini3b0417.pl -I all_PRO_${i}mM.top  -T $theta -O $PWtop

done

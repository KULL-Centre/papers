#!/bin/bash
#

export CUDA_HOME=/usr/local/cuda
export AMBERHOME=/temp1/Software/amber16
source /temp1/Software/amber16/amber.sh

pmemd.cuda -O -i min.mdin -p EXG.explicit.parm7 -c EXG.rst7 -o EXG.min.mdout -r EXG.min.rst7 -ref EXG.rst7 -cpin EXG.cpin


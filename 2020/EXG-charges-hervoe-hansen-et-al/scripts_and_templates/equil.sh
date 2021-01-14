#!/bin/bash
#

export CUDA_HOME=/usr/local/cuda
export AMBERHOME=/temp1/Software/amber16
source /temp1/Software/amber16/amber.sh

pmemd.cuda -O -i equil.mdin -p EXG.explicit.parm7 -c EXG.heat.rst7 -cpin EXG.cpin -o EXG.equil.mdout -r EXG.equil.rst7 -x EXG.equil.nc -cprestrt EXG.equil.cpin -cpout EXG.equil.cpout

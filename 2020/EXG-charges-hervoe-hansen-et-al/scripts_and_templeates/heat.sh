#!bin/bash
#

export CUDA_HOME=/usr/local/cuda
export AMBERHOME=/temp1/Software/amber16
source /temp1/Software/amber16/amber.sh

pmemd.cuda -O -i heat.mdin -c EXG.min.rst7 -p EXG.explicit.parm7 -cpin EXG.cpin -o EXG.heat.mdout -r EXGwt.heat.rst7 -x EXGwt.heat.nc


import os
import subprocess
import time
from jinja2 import Template
from analyse import *

proteins = initProteins()
proteins.to_pickle('proteins.pkl')

submission = Template("""#!/bin/bash
#SBATCH --job-name=m_{{name}}{{ff}}
#SBATCH --ntasks=1
#SBATCH --partition=qgpu
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=20000
#SBATCH -t 15:00:00
#SBATCH -o maps_{{name}}_{{ff}}_{{run}}.out
#SBATCH -e maps_{{name}}_{{ff}}_{{run}}.err

source /home/gitesei/.bashrc
conda activate hoomd
module purge
module load cmake/3.9.4 gcc/6.5.0 openmpi/4.0.3 llvm/7.0.0 cuda/9.2.148

start=$(date +%s.%N)

python ./maps_calc.py --name {{name}} --ff {{ff}} --run {{run}}

duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`

echo $execution_time""")

for name in ['ht40']:
    for i in [1,5,6]:
        for k in [1,2]:
            with open('maps_{:s}_ff{:d}0_{:d}.sh'.format(name,i,k), 'w') as submit:
                submit.write(submission.render(name=name,ff='{:d}0'.format(i),run=k))
            subprocess.run(['sbatch','--partition','qgpu','maps_{:s}_ff{:d}0_{:d}.sh'.format(name,i,k)])
            time.sleep(2)

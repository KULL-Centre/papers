import os
import subprocess
import time
from jinja2 import Template
from analyse import *

submission = Template("""#!/bin/bash
#SBATCH --job-name=PRE{{name}}{{ff}}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000
#SBATCH -t 15:00:00
#SBATCH -o maps_{{name}}_{{ff}}_{{run}}.out
#SBATCH -e maps_{{name}}_{{ff}}_{{run}}.err
#SBATCH --partition=sbinlab_ib

source /groups/sbinlab/giulio/.bashrc
conda activate hoomd

start=$(date +%s.%N)

python ./calc_map.py --name {{name}} --ff {{ff}} --run {{run}}

duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`

echo $execution_time""")

for name in ['A2','FUS']:
    for i in [1]:
        for k in [1]:
            with open('maps_{:s}_ff{:d}0_{:d}.sh'.format(name,i,k), 'w') as submit:
                submit.write(submission.render(name=name,ff='{:d}0'.format(i),run=k))
            subprocess.run(['sbatch','--partition','sbinlab','maps_{:s}_ff{:d}0_{:d}.sh'.format(name,i,k)])
            time.sleep(2)

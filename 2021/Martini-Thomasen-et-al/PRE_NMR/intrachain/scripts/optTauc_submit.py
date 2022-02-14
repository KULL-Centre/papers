import os
import subprocess
import time
from jinja2 import Template
from analyse import *

proteins = initProteins()
proteins.to_pickle('proteins.pkl')

submission = Template("""#!/bin/bash
#SBATCH --job-name=PRE{{name}}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH -t 20:00:00
#SBATCH -o PRE_{{name}}.out
#SBATCH -e PRE_{{name}}.err
#SBATCH --partition=sbinlab

source /groups/sbinlab/giulio/.bashrc
conda activate hoomd

start=$(date +%s.%N)

python ./optTauc.py --name {{name}}

duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`

echo $execution_time""")

for name in ['A2','FUS','OPN','aSyn_PRE','htau40']:
    with open('PRE_{}.sh'.format(name), 'w') as submit:
        submit.write(submission.render(name=name))
    subprocess.run(['sbatch','PRE_{}.sh'.format(name)])
    time.sleep(2)

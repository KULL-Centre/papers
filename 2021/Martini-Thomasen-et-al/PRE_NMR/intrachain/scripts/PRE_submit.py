import os
import subprocess
import time
from jinja2 import Template
from analyse import *

proteins = initProteins()
proteins.to_pickle('proteins.pkl')

submission = Template("""#!/bin/bash
#SBATCH --job-name=PRE{{name}}{{ff}}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=2000
#SBATCH -t 20:00:00
#SBATCH -o PRE_{{name}}_{{ff}}.out
#SBATCH -e PRE_{{name}}_{{ff}}.err
#SBATCH --partition=sbinlab_ib

source /groups/sbinlab/giulio/.bashrc
conda activate hoomd

start=$(date +%s.%N)

python ./calcPREs.py --name {{name}} --ff {{ff}}

duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`

echo $execution_time""")

lambda_range = [1.00,1.06,1.08]
for name in ['aSyn','FUS','A2']:
    for ff in lambda_range:
        with open('PRE_{:s}_{:.2f}.sh'.format(name,ff), 'w') as submit:
            submit.write(submission.render(name=name,ff='{:.2f}'.format(ff)))
        subprocess.run(['sbatch','--partition','sbinlab_ib','PRE_{:s}_{:.2f}.sh'.format(name,ff)])
        time.sleep(20)

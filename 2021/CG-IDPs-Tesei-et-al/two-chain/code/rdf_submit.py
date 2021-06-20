import os
import subprocess
import itertools
from jinja2 import Template
import time

submission = Template("""#!/bin/bash
#SBATCH --job-name=rdf-{{ff}}_{{name}}
#SBATCH --partition=sbinlab
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#SBATCH -t 2:00:00
#SBATCH -o rdf-{{ff}}_{{name}}.out
#SBATCH -e rdf-{{ff}}_{{name}}.err

source /groups/sbinlab/giulio/.bashrc
conda activate hoomd
python ./rdf_calc.py --name {{name}} --ff {{ff}}""")

for name,ff in itertools.product(['A2','FUS'],['M1','M2','M3']):
    with open('rdf-{:s}_{:s}.sh'.format(name,ff), 'w') as submit:
        submit.write(submission.render(name=name,ff=ff))
    subprocess.run(['sbatch','rdf-{:s}_{:s}.sh'.format(name,ff)])
    time.sleep(10)

import os
import subprocess
from jinja2 import Template

submission = Template("""#!/bin/bash
#SBATCH --job-name={{name}}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 20:00:00
#SBATCH -o calc.out
#SBATCH -e calc.err
#SBATCH --mem-per-cpu=20000
#SBATCH --partition=sbinlab_ib

source /groups/sbinlab/giulio/.bashrc
conda activate hoomd

python ./calc.py --model {{model}}""")

model = 'M1' # available models: 'M1', 'M2', 'M3', 'AVG', 'HPSUrry' and 'HPS'

with open('calc.sh', 'w') as submit:
    submit.write(submission.render(name='analysis',model=model))
subprocess.run(['sbatch','calc.sh'])

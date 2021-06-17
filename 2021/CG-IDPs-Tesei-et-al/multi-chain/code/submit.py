from analyse import *
import os
import subprocess
from jinja2 import Template

proteins = initProteins()
proteins.to_pickle('proteins.pkl')

submission = Template("""#!/bin/bash
#SBATCH --job-name={{name}}_{{temp}}
#SBATCH --nodes=1           
#SBATCH --partition=qgpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH -t 16:00:00
#SBATCH -o {{name}}_{{temp}}.out
#SBATCH -e {{name}}_{{temp}}.err

source /home/gitesei/.bashrc
conda activate hoomd
module purge
module load cmake/3.9.4 gcc/6.5.0 openmpi/4.0.3 llvm/7.0.0 cuda/9.2.148

echo $SLURM_CPUS_PER_TASK
echo $SLURM_JOB_NODELIST

python ./simulate.py --name {{name}} --temp {{temp}}""")

r = pd.read_pickle('residues.pkl')
r.lambdas = r['M1']
r.to_pickle('residues.pkl')

for name,prot in proteins.loc[['P4D','P7R','M6R','P7FM7Y','M8FP4Y','M10R']].iterrows():
    if not os.path.isdir(name):
        os.mkdir(name)
    for temp in [323]:
        if not os.path.isdir(name+'/{:d}'.format(temp)):
            os.mkdir(name+'/{:d}'.format(temp))
        with open('{:s}_{:d}.sh'.format(name,temp), 'w') as submit:
            submit.write(submission.render(name=name,temp='{:d}'.format(temp)))
        subprocess.run(['sbatch','{:s}_{:d}.sh'.format(name,temp)])

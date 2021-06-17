from analyse import *
import os
import subprocess
from jinja2 import Template

proteins = initProteinsDimers()
proteins.to_pickle('proteins.pkl')

submission = Template("""#!/bin/bash
#SBATCH --job-name={{name}}_{{ff}}
#SBATCH --nodes=1           
#SBATCH --partition=qgpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH -t 16:00:00
#SBATCH --mem-per-cpu=1500
#SBATCH -o {{name}}_{{ff}}_{{run}}.out
#SBATCH -e {{name}}_{{ff}}_{{run}}.err

source /home/gitesei/.bashrc
conda activate hoomd
module purge
module load cmake/3.9.4 gcc/6.5.0 openmpi/4.0.3 llvm/7.0.0 cuda/9.2.148

echo $SLURM_CPUS_PER_TASK
echo $SLURM_JOB_NODELIST

python ./simulate.py --name {{name}} --ff {{ff}} --run {{run}}""")

for name in ['ht40']:
    if not os.path.isdir(name):
        os.mkdir(name)
    for i in [1,2,3]:
        for j in range(1):
            if not os.path.isdir('{:s}/M{:d}{:d}'.format(name,i,j)):
                os.mkdir('{:s}/M{:d}{:d}'.format(name,i,j))
            for k in range(1,3):
                if not os.path.isdir('{:s}/M{:d}{:d}/run{:d}'.format(name,i,j,k)):
                    os.mkdir('{:s}/M{:d}{:d}/run{:d}'.format(name,i,j,k))
                with open('{:s}_M{:d}{:d}{:d}.sh'.format(name,i,j,k), 'w') as submit:
                    submit.write(submission.render(name=name,ff='{:d}{:d}'.format(i,j),run=k))
                subprocess.run(['sbatch','--partition','qgpu','{:s}_M{:d}{:d}{:d}.sh'.format(name,i,j,k)])

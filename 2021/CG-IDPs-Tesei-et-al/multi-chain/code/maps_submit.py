from analyse import *
import os
import subprocess
import shutil
from jinja2 import Template

proteins = initProteins()
proteins.to_pickle('proteins.pkl')

submission = Template("""#!/bin/bash
#SBATCH --job-name=m_{{name}}_{{temp}}
#SBATCH --nodes=1           
#SBATCH --partition=qgpu
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10000
#SBATCH --ntasks=1
#SBATCH -t 12:00:00
#SBATCH -o {{name}}_{{temp}}_maps.out
#SBATCH -e {{name}}_{{temp}}_maps.err

source /home/gitesei/.bashrc
conda activate hoomd
module purge
module load cmake/3.9.4 gcc/6.5.0 openmpi/4.0.3 llvm/7.0.0 cuda/9.2.148

python ./maps_calc.py --name {{name}} --temp {{temp}}""")

r = pd.read_pickle('residues.pkl')
r.lambdas = r['ff10']
r.to_pickle('residues.pkl')

for name,prot in proteins.loc[['A2','FUS','M4D','P4D','A1','P8D','A1NLS']].iterrows():
    for temp in [323]:
        with open('{:s}_{:d}_maps.sh'.format(name,temp), 'w') as submit:
            submit.write(submission.render(name=name,temp='{:d}'.format(temp)))
        subprocess.run(['sbatch','{:s}_{:d}_maps.sh'.format(name,temp)])

from analyse import *
import os
import subprocess
from jinja2 import Template
import time

submission_1 = Template("""#!/bin/sh
#SBATCH --job-name={{name}}_{{replica}}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=sbinlab_ib2
#SBATCH --mem=1GB
#SBATCH -t 10:00:00
#SBATCH -o {{name}}/{{cycle}}/err_{{replica}}
#SBATCH -e {{name}}/{{cycle}}/out_{{replica}}

source /groups/sbinlab/giulio/.bashrc

conda activate hoomd

echo $SLURM_CPUS_PER_TASK

echo $SLURM_CPUS_ON_NODE

python ./simulate.py --name {{name}} --cycle {{cycle}} --replica {{replica}} --cutoff {{cutoff}}""")

submission_2 = Template("""#!/bin/sh
#SBATCH --job-name={{name}}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=sbinlab_ib2
#SBATCH --mem=10GB
#SBATCH -t 1:00:00
#SBATCH -o {{name}}/{{cycle}}/merge_out
#SBATCH -e {{name}}/{{cycle}}/merge_err

source /groups/sbinlab/giulio/.bashrc

conda activate hoomd

echo $SLURM_CPUS_PER_TASK

echo $SLURM_CPUS_ON_NODE

python ./merge_replicas.py --name {{name}}""")

submission_3 = Template("""#!/bin/sh
#SBATCH --job-name=opt_M1
#SBATCH --nodes=1
#SBATCH --partition=sbinlab_ib2
#SBATCH --exclusive
#SBATCH --mem=90GB
#SBATCH -t 10:00:00
#SBATCH --dependency=afterok{% for id in jobid %}:{{id}}{% endfor %}
#SBATCH -o {{cycle}}_out
#SBATCH -e {{cycle}}_err

source /groups/sbinlab/giulio/.bashrc

conda activate hoomd

echo $SLURM_CPUS_PER_TASK

echo $SLURM_CPUS_ON_NODE

declare -a proteinsPRE_list=({{proteins}})

for name in ${proteinsPRE_list[@]}
do
cp -r ../expPREs/$name/expPREs $name
python ./pulchra.py --name $name --num_cpus $SLURM_CPUS_ON_NODE --pulchra /sbinlab/giulio/pulchra_306/pulchra
done

python ./optimize.py --log $SLURM_SUBMIT_DIR --cycle {{cycle}} --num_cpus $SLURM_CPUS_ON_NODE --cutoff {{cutoff}}""")

cycle = 12
cutoff = 2.0

proteins = initProteins(cycle)
proteinsRgs = initProteinsRgs(cycle)
allproteins = pd.concat((proteins,proteinsRgs),sort=True)
allproteins['N'] = allproteins['fasta'].apply(lambda x : len(x))
allproteins = allproteins.sort_values('N')

proteins.to_pickle('proteins.pkl')
proteinsRgs.to_pickle('proteinsRgs.pkl')
allproteins.to_pickle('allproteins.pkl')

prot_to_simulate = []
for name, prot in allproteins.iterrows():
    if not os.path.isfile(prot.path+'/{:s}.dcd'.format(name)):
        prot_to_simulate.append(name)
    else:
        n_frames = md.load(prot.path+'/{:s}.dcd'.format(name),
                           top=prot.path+'/{:s}.pdb'.format(name)).n_frames
        if n_frames!=5000:
            print(name,n_frames)
            prot_to_simulate.append(name)

print('Simulating {:d} sequences'.format(len(prot_to_simulate)))
print(prot_to_simulate)
jobid_1 = []
"""
for name in prot_to_simulate:
    if not os.path.isdir(name):
        os.mkdir(name)
    if not os.path.isdir(name+'/{:d}'.format(cycle)):
        os.mkdir(name+'/{:d}'.format(cycle))
    for replica in range(10):
        if not os.path.isfile('{:s}/{:d}/{:d}.gsd'.format(name,cycle,replica)):
            with open('{:s}_{:d}_{:d}.sh'.format(name,cycle,replica), 'w') as submit:
                submit.write(submission_1.render(name=name,cycle='{:d}'.format(cycle),
                                                 replica='{:d}'.format(replica),cutoff='{:.1f}'.format(cutoff)))
            proc = subprocess.run(['sbatch','{:s}_{:d}_{:d}.sh'.format(name,cycle,replica)],capture_output=True)
            jobid_1.append(int(proc.stdout.split(b' ')[-1].split(b'\\')[0]))
        elif os.path.getsize('{:s}/{:d}/{:d}.gsd'.format(name,cycle,replica)) < 400000:
            with open('{:s}_{:d}_{:d}.sh'.format(name,cycle,replica), 'w') as submit:
                submit.write(submission_1.render(name=name,cycle='{:d}'.format(cycle),
                                                 replica='{:d}'.format(replica),cutoff='{:.1f}'.format(cutoff)))
            proc = subprocess.run(['sbatch','{:s}_{:d}_{:d}.sh'.format(name,cycle,replica)],capture_output=True)
            jobid_1.append(int(proc.stdout.split(b' ')[-1].split(b'\\')[0]))
        time.sleep(0.5)
"""
jobid_2 = []
for name in prot_to_simulate:
    if not os.path.isdir(name):
        os.mkdir(name)
    if not os.path.isdir(name+'/{:d}'.format(cycle)):
        os.mkdir(name+'/{:d}'.format(cycle))
    with open('{:s}_{:d}.sh'.format(name,cycle), 'w') as submit:
        submit.write(submission_2.render(jobid=jobid_1,name=name,cycle='{:d}'.format(cycle)))
    proc = subprocess.run(['sbatch','{:s}_{:d}.sh'.format(name,cycle)],capture_output=True)
    jobid_2.append(int(proc.stdout.split(b' ')[-1].split(b'\\')[0]))
    time.sleep(0.5)
with open('opt_{:d}.sh'.format(cycle), 'w') as submit:
    submit.write(submission_3.render(jobid=jobid_2,proteins=' '.join(proteins.index),
                              cycle='{:d}'.format(cycle),cutoff='{:.1f}'.format(cutoff)))
subprocess.run(['sbatch','opt_{:d}.sh'.format(cycle)])

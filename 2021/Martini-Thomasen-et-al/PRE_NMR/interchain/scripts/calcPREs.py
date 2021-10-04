from analyse import *
import MDAnalysis
import time
import os
import glob
import sys
from DEERPREdict.PRE import PREpredict
from DEERPREdict.utils import Operations
from argparse import ArgumentParser
import ray
import logging
import psutil
import shutil

parser = ArgumentParser()
parser.add_argument('--name',dest='name',type=str,required=True)
parser.add_argument('--ff',dest='ff',type=str,required=True)
parser.add_argument('--run',dest='run',type=int,required=True)
args = parser.parse_args()

proteins = pd.read_pickle('proteins.pkl')

# create directory called calcPREs 
# the calculated PREs will be saved in this directory
for _, prot in proteins.loc[[args.name]].iterrows():
    if not os.path.isdir(prot.path+'{:s}/run{:d}/calcPREs'.format(args.ff,args.run)):
        os.mkdir(prot.path+'{:s}/run{:d}/calcPREs'.format(args.ff,args.run))

# determine number of label positions for each protein
proc_PRE = [(label,name) for name,prot in proteins.loc[[args.name]].iterrows() for label in prot.labels]
# use total number of label positions as number of cpus
num_cpus = len(proc_PRE)
     
ray.init(num_cpus=num_cpus)

@ray.remote
def evaluatePRE(n, label, name, prot):
    prefix = prot.path+'{:s}/run{:d}/calcPREs/res'.format(args.ff,args.run)
    filename = prefix+'-{:d}.pkl'.format(label)
    u = MDAnalysis.Universe(prot.path+'{:s}/run{:d}/allatom.pdb'.format(args.ff,args.run),
            prot.path+'{:s}/run{:d}/allatom.dcd'.format(args.ff,args.run))
    load_file = False
    PRE = PREpredict(u, label, chains=['A','B'], log_file = prot.path+'{:s}/run{:d}/log'.format(args.ff,args.run),
                     temperature = prot.temp, atom_selection = 'N', sigma_scaling = 1.0)
    PRE.run(output_prefix = prefix+'AB', load_file = load_file, tau_t = 1e-10, tau_c = 1e-09, r_2 = 10, wh = prot.wh)
    PRE = PREpredict(u, label, chains=['B','A'], log_file = prot.path+'{:s}/run{:d}/log'.format(args.ff,args.run),
                     temperature = prot.temp, atom_selection = 'N', sigma_scaling = 1.0)
    PRE.run(output_prefix = prefix+'BA', load_file = load_file, tau_t = 1e-10, tau_c = 1e-09, r_2 = 10, wh = prot.wh)

for name,prot in proteins.loc[[args.name]].iterrows():
   s = md.load(prot.path+'{:s}/run{:d}/allatom.gro'.format(args.ff,args.run))
   a = s.atom_slice(s.top.select('index 0 to {:d}'.format(int(s.n_atoms/2-1))))
   b = s.atom_slice(s.top.select('index {:d} to {:d}'.format(int(s.n_atoms/2),s.n_atoms)))
   stacked_s = a.stack(b)
   stacked_s.save_pdb(prot.path+'{:s}/run{:d}/allatom.pdb'.format(args.ff,args.run))
   traj = md.load(prot.path+'{:s}/run{:d}/allatom.xtc'.format(args.ff,args.run),
           top = prot.path+'{:s}/run{:d}/allatom.pdb'.format(args.ff,args.run))
   traj = traj.image_molecules(inplace=False, anchor_molecules=[set(traj.top.chain(0).atoms)],
           other_molecules=[set(traj.top.chain(1).atoms)], make_whole=False)
   traj[0].save_pdb(prot.path+'{:s}/run{:d}/allatom.pdb'.format(args.ff,args.run))
   traj.save_dcd(prot.path+'{:s}/run{:d}/allatom.dcd'.format(args.ff,args.run))

# run PRE calculations
time0 = time.time()
ray.get([evaluatePRE.remote(n,label,name,proteins.loc[name]) for n,(label,name) in enumerate(proc_PRE)])
print(args.name,args.ff,args.run)
print('Timing {:.3f}'.format(time.time()-time0))

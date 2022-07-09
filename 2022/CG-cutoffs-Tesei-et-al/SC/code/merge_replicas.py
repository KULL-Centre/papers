from analyse import *
import hoomd
import hoomd.md
from hoomd import azplugins
import time
import itertools
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--name',nargs='?',const='', type=str)
args = parser.parse_args()

print(hoomd.__file__)

def centerDCD(residues,name,prot):
    top = md.Topology()
    chain = top.add_chain()
    for resname in prot.fasta:
        residue = top.add_residue(residues.loc[resname,'three'], chain)
        top.add_atom(residues.loc[resname,'three'], element=md.element.carbon, residue=residue)
    for i in range(len(prot.fasta)-1):
        top.add_bond(top.atom(i),top.atom(i+1))
    traj = md.load(prot.path+'/0.gsd')[100:]
    traj.top = top
    for i in range(1,10):
        t = md.load(prot.path+'/{:d}.gsd'.format(i))[100:]
        t.top = top
        traj = md.join([traj,t])
    traj = traj.image_molecules(inplace=False, anchor_molecules=[set(traj.top.chain(0).atoms)], make_whole=True)
    traj.center_coordinates()
    traj.xyz += traj.unitcell_lengths[0,0]/2
    print('Number of frames: {:d}'.format(traj.n_frames))
    traj.save_dcd(prot.path+'/{:s}.dcd'.format(name))
    traj[0].save_pdb(prot.path+'/{:s}.pdb'.format(name))
    for i in range(10):
        os.remove(prot.path+'/{:d}.gsd'.format(i))

residues = pd.read_csv('residues.csv').set_index('one',drop=False)

allproteins = pd.read_pickle('allproteins.pkl')

t0 = time.time()
centerDCD(residues,args.name,allproteins.loc[args.name])
print('Timing Simulation {:.3f}'.format(time.time()-t0))

import numpy as np
import pandas as pd
import itertools
import mdtraj as md
from mdtraj import element
import pickle
import time
import os
from argparse import ArgumentParser

def individual_rdfs(protname,run,ff,width,proteins,residues):
    """
    Calculate radial distribution function between the 
    centers of mass of two protein chains from a single traj. 
    width is the bin width
    of the histogram used to construct the graph.
    """
    
    upper = proteins.loc[protname].L / 2 # upper limit of the range of distances (in nm)
    n_chains = 2
    masses = residues.loc[proteins.loc[protname].fasta,'MW'].values
    masses[0] += 2
    masses[-1] += 16
    radii = residues.loc[proteins.loc[protname].fasta,'sigmas'].values/2
    # define topology that includes bead masses to calculate correct center of mass
    top = md.Topology()
    for _ in range(n_chains):
        chain = top.add_chain()
        residue = top.add_residue('C{:d}'.format(chain.index), chain, resSeq=chain.index)
        for i,resname in enumerate(proteins.loc[protname].fasta):
            element.Element._elements_by_symbol.pop('A'+resname, None)
            el = element.Element.__new__(element.Element, 1, 'A'+resname, 'A'+resname, masses[i], radii[i])
            atom = top.add_atom('A'+resname, element=el, residue=residue)
        for i in range(chain.n_atoms-1):
            top.add_bond(chain.atom(i),chain.atom(i+1))
   
    # load trajectory data 
    t = md.load('{:s}/{:s}/run{:d}/{:s}.dcd'.format(protname,ff,run,protname),top=top)
    t = t[66:]

    # create trajectory and topology for centers of mass
    cmtop = md.Topology()
    cmpos = []
    for chain in t.top.chains:
        chain = cmtop.add_chain()
        res = cmtop.add_residue('CM', chain, resSeq=chain.index)
        cmtop.add_atom('CM', element=t.top.atom(0).element, residue=res)
        cmpos.append(md.compute_center_of_mass(
            t.atom_slice(t.top.select('chainid {:d}'.format(chain.index)))))
    cmpos = np.swapaxes(np.array(cmpos),0,1)
    cmtraj = md.Trajectory(cmpos, cmtop, t.time, t.unitcell_lengths, t.unitcell_angles)
    # calculate the rdf between the centers of mass
    for k in range(4): # divide each replica into 4 chunks
        b = 29150*k
        e = b+29150
        r,rdf = md.compute_rdf(cmtraj[b:e], [[0,1]], r_range=(0,upper), bin_width = width, periodic=True)
        # save results
        np.savetxt('rdfs/{:s}_{:s}_{:d}_{:d}.dat'.format(protname,ff,run,k),np.c_[r,rdf])
    
def concatenated_rdf(protname,n_runs,ff,width,proteins,residues):
    """
    Caculate rdf from a single, long trajectory after concatenation 
    of several trajectories. The concatenation takes place after
    the center-of-mass trajectories have been constructed.
    The calculation will be performed for the first n_runs trajectories. 
    """
    
    # define topology that includes bead masses to calculate correct center of mass
    upper = proteins.loc[protname].L / 2
    n_chains = 2
    masses = residues.loc[proteins.loc[protname].fasta,'MW'].values
    masses[0] += 2
    masses[-1] += 16
    radii = residues.loc[proteins.loc[protname].fasta,'sigmas'].values/2
    top = md.Topology()
    for _ in range(n_chains):
        chain = top.add_chain()
        residue = top.add_residue('C{:d}'.format(chain.index), chain, resSeq=chain.index)
        for i,resname in enumerate(proteins.loc[protname].fasta):
            element.Element._elements_by_symbol.pop('A'+resname, None)
            el = element.Element.__new__(element.Element, 1, 'A'+resname, 'A'+resname, masses[i], radii[i])
            atom = top.add_atom('A'+resname, element=el, residue=residue)
        for i in range(chain.n_atoms-1):
            top.add_bond(chain.atom(i),chain.atom(i+1))

    cmtrajs = []    
    for run in range(1,n_runs+1):
        # load trajectory data 
        t = md.load('{:s}/{:s}/run{:d}/{:s}.dcd'.format(protname,ff,run,protname),top=top)
        t = t[66:]
        # create trajectory and topology for centers of mass
        cmtop = md.Topology()
        cmpos = []
        for chain in t.top.chains:
            chain = cmtop.add_chain()
            res = cmtop.add_residue('CM', chain, resSeq=chain.index)
            cmtop.add_atom('CM', element=t.top.atom(0).element, residue=res)
            cmpos.append(md.compute_center_of_mass(
                t.atom_slice(t.top.select('chainid {:d}'.format(chain.index)))))
        cmpos = np.swapaxes(np.array(cmpos),0,1)
        cmtraj = md.Trajectory(cmpos, cmtop, t.time, t.unitcell_lengths, t.unitcell_angles)
        cmtrajs.append(cmtraj)
    # create pseudo-traj based on all individual trajectories 
    pseudo_traj = md.join(cmtrajs)
    # calculate the rdf between the centers of mass
    r,rdf = md.compute_rdf(pseudo_traj, [[0,1]], r_range=(0,upper), bin_width = width, periodic=True)
    # save results
    np.savetxt('rdfs/{:s}_{:s}.dat'.format(protname,ff),np.c_[r,rdf])

parser = ArgumentParser()
parser.add_argument('--name',dest='name',type=str,required=True)
parser.add_argument('--ff',dest='ff',type=str,required=True)
args = parser.parse_args()    

proteins = pd.read_pickle('proteins.pkl')
residues = pd.read_pickle('residues.pkl')
residues = residues.set_index('one')

# create directory to save files in
if not os.path.isdir('rdfs'):
    os.mkdir('rdfs')

# calculate individual rdfs for blocks of 875 ns
for run in range(1,11):
    individual_rdfs(args.name,run,args.ff,0.1,proteins,residues)

# calculate rdf using data from all trajs
concatenated_rdf(args.name,10,args.ff,0.1,proteins,residues)

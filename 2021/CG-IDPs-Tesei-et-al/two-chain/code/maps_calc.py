from analyse import *
import MDAnalysis
import time
import os
import glob
import sys
import pandas as pd
import numpy as np
import mdtraj as md
import itertools
from mdtraj import element
from argparse import ArgumentParser
from scipy.optimize import least_squares
import time

parser = ArgumentParser()
parser.add_argument('--name',nargs='?',const='', type=str)
parser.add_argument('--ff',nargs='?',const='', type=str)
parser.add_argument('--run',nargs='?',const='', type=int)
args = parser.parse_args()

def calcContactMap(traj,df,prot,name,contact_map,energies,lambdamap,sigmamap,lj_eps):
    sel1 = traj.top.select('chainid 0')
    sel2 = traj.top.select('chainid 1')
    pairs_indices = traj.top.select_pairs(sel1,sel2)
    d = md.compute_distances(traj,pairs_indices)

    d[d>4] = np.inf
    lj = lambda x,sig : 4*lj_eps*((sig/x)**12-(sig/x)**6)
    ah = lambda x,sig,l : np.where(x<=np.power(2.,1./6)*sig,lj(x,sig)+lj_eps*(1-l),l*lj(x,sig))
    en_ah = np.apply_along_axis(lambda a: ah(a,sigmamap,lambdamap), 1, d)
    
    energies = np.append(energies, en_ah.sum(axis=1))
    contact_map += en_ah.mean(axis=0)
    print(energies.shape)
    print(d.shape)
    del d
    del en_ah
    return contact_map, energies

def trajCM(df,proteins,name,ff,run):
    # this function finds the index of the chain at the center of the slab for each frame
    path = '{:s}/ff{:s}/run{:d}/'.format(name,ff,run)
    prot = proteins.loc[name]
    masses = df.loc[prot.fasta,'MW'].values
    masses[0] += 2                                                                                                   
    masses[-1] += 16
    radii = df.loc[prot.fasta,'sigmas'].values/2

    t = md.load_pdb(path+'{:s}.pdb'.format(name))

    n_chains = int( t.n_atoms / len(prot.fasta) )

    n_chains = 2
    top = md.Topology()
    for _ in range(n_chains):
        chain = top.add_chain()
        residue = top.add_residue('C{:d}'.format(chain.index), chain, resSeq=chain.index)
        for i,resname in enumerate(proteins.loc[name].fasta):
            element.Element._elements_by_symbol.pop('A'+resname, None)
            el = element.Element.__new__(element.Element, 1, 'A'+resname, 'A'+resname, masses[i], radii[i])
            atom = top.add_atom('A'+resname, element=el, residue=residue)
        for i in range(chain.n_atoms-1):
            top.add_bond(chain.atom(i),chain.atom(i+1))
   
    # load trajectory data 
    t = md.load(path+'{:s}.dcd'.format(name),top=top)
    t = t[66:]

    pairs = np.array(list(itertools.product(prot.fasta,prot.fasta)))
    pairs = np.core.defchararray.add(pairs[:,0],pairs[:,1])
    _, lj_eps, lj_lambda, lj_sigma, _, _ = genParamsLJ(df,name,prot)
    dflambda = lj_lambda.unstack()
    dflambda.index = dflambda.index.map('{0[0]}{0[1]}'.format)
    dfsigma = lj_sigma.unstack()
    dfsigma.index = dfsigma.index.map('{0[0]}{0[1]}'.format)

    lambdamap = dflambda.loc[pairs]
    sigmamap = dfsigma.loc[pairs]
    print(t.n_frames,lambdamap.shape,sigmamap.shape,lj_eps)

    contact_map = np.zeros(len(prot.fasta)*len(prot.fasta))
    energies = np.empty(0)
    #for traj in :
    n_chunks = 8
    len_chunk = t.n_frames//n_chunks
    print(t.n_frames,len_chunk*n_chunks,lambdamap.shape,sigmamap.shape,lj_eps)
    for k in range(n_chunks):
        b = len_chunk*k
        e = b+len_chunk
        contact_map, energies = calcContactMap(t[b:e],df,prot,name,contact_map,energies,lambdamap,sigmamap,lj_eps)
        print('traj',t[b:e].n_frames)
    contact_map = contact_map.reshape(len(prot.fasta),len(prot.fasta)) / n_chunks

    return contact_map, energies

df = pd.read_pickle('residues.pkl').set_index('one')
df.lambdas = df['ff{:s}'.format(args.ff)]
proteins = pd.read_pickle('proteins.pkl')
prot = proteins.loc[args.name]

t0 = time.time()

contact_map, energies = trajCM(df,proteins,args.name,args.ff,args.run)
np.save('maps/{:s}_{:s}_{:d}_map.npy'.format(args.name,args.ff,args.run),contact_map)
np.save('maps/{:s}_{:s}_{:d}_energies.npy'.format(args.name,args.ff,args.run),energies)

print('Timing {:.3f}'.format(time.time()-t0))

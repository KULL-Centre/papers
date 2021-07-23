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
parser.add_argument('--temp',nargs='?',const='', type=int)
args = parser.parse_args()

lj = lambda x,sig,lj_eps : 4*lj_eps*((sig/x)**12-(sig/x)**6)
ah = lambda x,sig,l,lj_eps : np.where(x<=np.power(2.,1./6)*sig,lj(x,sig,lj_eps)+lj_eps*(1-l),l*lj(x,sig,lj_eps))
 
def calcContactMap(cmap,ecmap,traj,chain1,chain2,lambdamap,sigmamap,lj_eps):
    sel1 = traj.top.select('chainid {:d}'.format(chain1))
    sel2 = traj.top.select('chainid {:d}'.format(chain2))
    pairs_indices = traj.top.select_pairs(sel1,sel2)
    d = md.compute_distances(traj,pairs_indices)
    d[d>4] = np.inf
    en_ah = np.apply_along_axis(lambda a: ah(a,sigmamap,lambdamap,lj_eps), 1, d)
    cmap += np.where(d<0.8,1,0).sum(axis=0)
    ecmap += en_ah.sum(axis=0)
    del d
    del en_ah
    return cmap, ecmap

def calcWidth(path,name,temp):
    # this function finds the z-positions that delimit the slab and the dilute phase
    h = np.load('{:s}_{:d}.npy'.format(name,temp))
    lz = (h.shape[1]+1)
    edges = np.arange(-lz/2.,lz/2.,1)/10
    dz = (edges[1]-edges[0])/2.
    z = edges[:-1]+dz
    profile = lambda x,a,b,c,d : .5*(a+b)+.5*(b-a)*np.tanh((np.abs(x)-c)/d)
    residuals = lambda params,*args : ( args[1] - profile(args[0], *params) )
    hm = np.mean(h,axis=0)
    z1 = z[z>0]
    h1 = hm[z>0]
    z2 = z[z<0]
    h2 = hm[z<0]
    p0=[hm.min(),hm.max(),3,1]
    res1 = least_squares(residuals, x0=p0, args=[z1, h1], bounds=([0]*4,[1e3]*4))
    res2 = least_squares(residuals, x0=p0, args=[z2, h2], bounds=([0]*4,[1e3]*4))
    cutoff1 = .5*(np.abs(res1.x[2]-.5*res1.x[3])+np.abs(-res2.x[2]+.5*res2.x[3]))
    cutoff2 = .5*(np.abs(res1.x[2]+6*res1.x[3])+np.abs(-res2.x[2]-6*res2.x[3]))
    return cutoff1, cutoff2

def trajCM(df,proteins,name,temp):
    # this function finds the index of the chain at the center of the slab for each frame
    path = '{:s}/{:d}/'.format(name,temp)
    cutoff1, cutoff2 = calcWidth(path,name,temp)
    print(name,cutoff1,cutoff2)
    prot = proteins.loc[name]
    masses = df.loc[prot.fasta,'MW'].values
    masses[0] += 2                                                                                                   
    masses[-1] += 16
    radii = df.loc[prot.fasta,'sigmas'].values/2

    t = md.load_pdb(path+'top.pdb')

    n_chains = int( t.n_atoms / len(prot.fasta) )

    top = md.Topology()
    for _ in range(n_chains):
        chain = top.add_chain()
        residue = top.add_residue('C{:d}'.format(chain.index), chain, resSeq=chain.index)
        for i,resname in enumerate(prot.fasta):
            # add an element with unique name to the dictionary. the letter A is prepended to avoid doubles (e.g. cysteine and carbon)
            element.Element._elements_by_symbol.pop('A'+resname, None)
            el = element.Element.__new__(element.Element, 1, 'A'+resname, 'A'+resname, masses[i], radii[i])
            atom = top.add_atom('A'+resname, element=el, residue=residue)
        for i in range(chain.n_atoms-1):
            top.add_bond(chain.atom(i),chain.atom(i+1))

    t = md.load_dcd(path+'traj.dcd',top)
    t.xyz -= t.unitcell_lengths[0,:]/2
    t.make_molecules_whole(inplace=True)
    t = t[1200:] # skip first 0.3 us

    print(t.n_frames)

    rgsC = np.empty(0) # Rgs of chains in the slab
    rgsD = np.empty(0) # Rgs of chains in dilute phase
    indices = np.zeros((n_chains,t.n_frames))
    middle_dist = np.zeros((n_chains,t.n_frames))
    for res in t.top.residues:
        #print('chain',res.index)
        chainslice = t.atom_slice(t.top.select('chainid {:d}'.format(res.index)))
        chaincm = md.compute_center_of_mass(chainslice)
        mask_in = np.abs(chaincm[:,2])<cutoff1
        mask_out = np.abs(chaincm[:,2])>cutoff2
        rgsC = np.append(rgsC, md.compute_rg(chainslice,masses=masses)[mask_in])
        rgsD = np.append(rgsD, md.compute_rg(chainslice,masses=masses)[mask_out])
        indices[res.index] = mask_in
        middle_dist[res.index] = np.abs(chaincm[:,2])

    middle_chain = np.argmin(middle_dist,axis=0) # indices of chains at the center of the slab
    del middle_dist

    pairs = np.array(list(itertools.product(prot.fasta,prot.fasta)))
    pairs = np.core.defchararray.add(pairs[:,0],pairs[:,1])
    _, lj_eps, lj_lambda, lj_sigma, _, _, _ = genParamsLJ(df,name,prot)
    dflambda = lj_lambda.unstack()
    dflambda.index = dflambda.index.map('{0[0]}{0[1]}'.format)
    dfsigma = lj_sigma.unstack()
    dfsigma.index = dfsigma.index.map('{0[0]}{0[1]}'.format)
    lambdamap = dflambda.loc[pairs]
    sigmamap = dfsigma.loc[pairs]

    cmap = np.zeros(len(prot.fasta)*len(prot.fasta))
    ecmap = np.zeros(len(prot.fasta)*len(prot.fasta))
    for chain1 in np.unique(middle_chain):
        for chain2 in np.setdiff1d(np.arange(n_chains),[chain1]):
            ndx = ((middle_chain==chain1)*indices[chain2]).astype(bool)
            if np.any(ndx):
                cmap, ecmap = calcContactMap(cmap,ecmap,t[ndx],chain1,chain2,lambdamap,sigmamap,lj_eps)
    cmap = cmap.reshape(len(prot.fasta),len(prot.fasta)) / t.n_frames
    ecmap = ecmap.reshape(len(prot.fasta),len(prot.fasta)) / t.n_frames
    return cmap, ecmap, rgsC.flatten(), rgsD.flatten()

df = pd.read_csv('residues.csv').set_index('three',drop=False).set_index('one')
proteins = pd.read_pickle('proteins.pkl')
prot = proteins.loc[args.name]

t0 = time.time()

cmap, ecmap, rgsC, rgsD = trajCM(df,proteins,args.name,args.temp)
np.savetxt('{:s}_{:d}_rgC.dat'.format(args.name,args.temp),rgsC)
np.savetxt('{:s}_{:d}_rgD.dat'.format(args.name,args.temp),rgsD)
np.save('{:s}_{:d}_cmap.npy'.format(args.name,args.temp),cmap)
np.save('{:s}_{:d}_ecmap.npy'.format(args.name,args.temp),ecmap)
print('Timing {:.3f}'.format(time.time()-t0))

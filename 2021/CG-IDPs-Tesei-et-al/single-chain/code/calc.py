import MDAnalysis
import time
import os
import glob
import sys
import psutil
import logging
import shutil
import h5py
import pandas as pd
import numpy as np
import mdtraj as md
import itertools
from mdtraj import element
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
import time
from argparse import ArgumentParser

proteins = pd.read_pickle('proteins.pkl')

def genParamsLJ(df,name,prot):
    fasta = prot.fasta.copy()
    r = df.copy()
    r.loc['X'] = r.loc[fasta[0]]
    r.loc['Z'] = r.loc[fasta[-1]]
    fasta[0] = 'X'
    fasta[-1] = 'Z'
    types = list(np.unique(fasta))
    MWs = [r.loc[a,'MW'] for a in types]
    sigmamap = pd.DataFrame((r.sigmas.values+r.sigmas.values.reshape(-1,1))/2,
                            index=r.sigmas.index,columns=r.sigmas.index)
    lambdamap = pd.DataFrame((r.lambdas.values+r.lambdas.values.reshape(-1,1))/2,
                             index=r.lambdas.index,columns=r.lambdas.index)
    lj_eps = 0.2*4.184
    pairs = np.array(list(itertools.combinations_with_replacement(types,2)))
    return pairs, lj_eps, lambdamap, sigmamap, fasta, types, MWs

def calcContactMap(d,df,prot,name):
    d[d>4] = np.inf
    pairs = np.array(list(itertools.combinations(prot.fasta,2)))
    pairs = np.core.defchararray.add(pairs[:,0],pairs[:,1])
    _, lj_eps, lj_lambda, lj_sigma, _, _, _ = genParamsLJ(df,name,prot)
    dflambda = lj_lambda.unstack()
    dflambda.index = dflambda.index.map('{0[0]}{0[1]}'.format)
    dfsigma = lj_sigma.unstack()
    dfsigma.index = dfsigma.index.map('{0[0]}{0[1]}'.format)

    lj = lambda x,sig : 4*lj_eps*((sig/x)**12-(sig/x)**6)
    ah = lambda x,sig,l : np.where(x<=np.power(2.,1./6)*sig,lj(x,sig)+lj_eps*(1-l),l*lj(x,sig))
    en_ah = np.apply_along_axis(lambda a: ah(a,dfsigma.loc[pairs],dflambda.loc[pairs]), 1, d)
    
    contact_map = en_ah.mean(axis=0)
    del d
    del en_ah
    del dfsigma
    del dflambda
    return contact_map

def calcRs(traj):
    pairs = traj.top.select_pairs('all','all')
    d = md.compute_distances(traj,pairs)
    dmap = np.where(d<0.8,1,0).mean(axis=0)
    dmean = d.mean(axis=0)
    ij = np.array(range(1,traj.n_atoms))
    diff = [x[1]-x[0] for x in pairs]
    dij = np.empty(0)
    for i in ij:
        dij = np.append(dij,dmean[diff==i].mean())
    return ij, dij, d, dmap 

def calcRg(df,traj):
    residues = [res.name for res in traj.top.atoms]
    masses = df.loc[residues,'MW'].values
    rgarray = md.compute_rg(traj,masses=masses)
    rg = np.sqrt( np.power(rgarray, 2).mean() )
    rgE = np.std(rgarray)/np.sqrt(rgarray.size)
    return rgarray, rg, rgE

parser = ArgumentParser()
parser.add_argument('--model',nargs='?',const='', type=str)
args = parser.parse_args()

residues = pd.read_csv('residues.csv',index_col='three')
residues.lambdas = residues[args.model]
t0 = time.time()

data = {}
f = lambda x,R0,v : R0*np.power(x,v)
for name in proteins.index:
    traj = md.load_dcd("{:s}/{:s}.dcd".format(name,name),"{:s}/{:s}.pdb".format(name,name))
    rgarray, _, _ = calcRg(residues,traj)
    traj = traj[5000:]
    ij, dij, d, dmap = calcRs(traj)
    popt, pcov = curve_fit(f,ij[ij>10],dij[ij>10],p0=[.4,.5])
    _, rg, rgE = calcRg(residues,traj)
    contact_map = calcContactMap(d,residues.set_index('one'),proteins.loc[name],name)
    data[name] = {'ij':ij, 'dij':dij, 'simNu':popt[1], 'simNuE':np.sqrt(np.diag(pcov))[1],
            'simRg':rg, 'simRgE':rgE, 'rgarray':rgarray, 'map':contact_map} 

df = pd.DataFrame(data=data).T

df.to_pickle('calc.pkl')
print('Timing {:.3f}'.format(time.time()-t0))

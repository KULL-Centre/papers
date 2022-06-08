import mdtraj as md
import numpy as np
import sys

def rh_kirk(n, conf):
    invrij = (1-1/n)*(1/md.compute_distances(conf,conf.top.select_pairs('all','all'))).mean(axis=1)
    return 1/invrij


name = sys.argv[1]
folder = '/storage1/francesco/PROJECTS/ENSEMBLES/TSCL-M1/'+name+'/CG_frames/'

rh = []
for i in range(0,20000):
    conf = md.load_pdb(folder+'frame{}.pdb'.format(i))
    n = len(list(conf.topology.residues))
    rh.append(rh_kirk(n,conf))

np.savetxt(name+'/Rh_Kirk_Ca.dat', rh)

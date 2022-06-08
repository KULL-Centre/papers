import mdtraj as md
import numpy as np
import sys

def rh_nygaard(n,rg):
    a1 = 0.216
    a2 = 4.06
    a3 = 0.821

    rg = rg*10

    num = a1 * (rg - a2*(n)**0.33)
    den = (n**0.6) - (n**0.33)
    rh_inv = ((num/den)+a3)/rg

    return(1/rh_inv)

name = sys.argv[1]
folder = '/storage1/francesco/PROJECTS/ENSEMBLES/TSCL-M1/'+name+'/CG_frames/'

rh = []
for i in range(0,20000):
    conf = md.load_pdb(folder+'frame{}.pdb'.format(i))
    n = len(list(conf.topology.residues))
    rg = md.compute_rg(conf)
    rh.append(rh_nygaard(n,rg))

np.savetxt(name+'/Rh_Nyg.dat', rh)

import mdtraj as md
import sys
import numpy as np

xtc = sys.argv[1]
pdb = sys.argv[2]

traj = md.load_xtc(xtc, top=pdb)
pairs = traj.top.select_pairs('all','all')
d = md.compute_distances(traj,pairs)
print(np.max(d))

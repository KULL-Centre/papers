import numpy as np
import pandas as pd
import mdtraj as md
import time
import ray
import os
import subprocess
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--num_cpus',dest='num_cpus',type=int)
parser.add_argument('--name',dest='name',type=str)
parser.add_argument('--pulchra',dest='pulchra_path',type=str)
args = parser.parse_args()

ray.init(num_cpus=args.num_cpus)

def fix_topology(dcd,pdb):
    """
    Changes atom names to CA
    """
    t = md.load_dcd(dcd,pdb)
    cgtop = md.Topology()
    cgchain = cgtop.add_chain()
    for atom in t.top.atoms:
        cgres = cgtop.add_residue(atom.name, cgchain)
        cgtop.add_atom('CA', element=md.element.carbon, residue=cgres)
    traj = md.Trajectory(t.xyz, cgtop, t.time, t.unitcell_lengths, t.unitcell_angles)
    traj = traj.superpose(traj, frame=0)
    return traj

@ray.remote
def run_pulchra(prot_path,pulchra_path,i,frame):
    """
    This function runs Pulchra on a single frame.
    """
    name = prot_path+'/{:d}.pdb'.format(i)
    frame.save(name)
    FNULL = open(os.devnull, 'w')
    subprocess.run([pulchra_path,name],stdout=FNULL,stderr=FNULL)    
    outname = prot_path+'/{:d}.rebuilt.pdb'.format(i)
    trajtemp = md.load(outname)
    os.remove(name)
    os.remove(outname)
    return trajtemp.xyz

def reconstruct_pulchra(pulchra_path,prot_name,prot,num_cpus):
    """
    This function reconstructs an all-atom trajectory from a Calpha trajectory.
    Input: trajectory nvt.xtc and topology nvt.gro file.
    n_procs: number of processors to use in parallel.
    Return: A reconstructed mdtraj trajectory.
    """
    t = fix_topology(prot.path+'/{:s}.dcd'.format(prot_name),prot.path+'/{:s}.pdb'.format(prot_name))
    name = prot.path+'/0.pdb'
    t[0].save(name)
    subprocess.run([pulchra_path,name])    
    s = md.load_pdb(prot.path+'/0.rebuilt.pdb')
    # n_blocks = t.n_frames // num_cpus
    xyz = np.empty((0,s.n_atoms,3))
    xyz = np.append( xyz, s.xyz )
    num_cpus = num_cpus - 1
    for j in range(1, t.n_frames, num_cpus):
        n = j+num_cpus if j+num_cpus<t.n_frames else t.n_frames
        xyz = np.append( xyz, np.vstack(ray.get([run_pulchra.remote(prot.path,pulchra_path,i,t[i]) for i in range(j,n)])) )
    allatom0 = md.Trajectory(xyz.reshape(t.n_frames,s.n_atoms,3), s.top, t.time, t.unitcell_lengths, t.unitcell_angles)
    top = md.Topology()
    chain = top.add_chain()
    for residue in allatom0.top.residues:
        res = top.add_residue(residue.name, chain, resSeq=residue.index+1)
        for atom in residue.atoms:
            top.add_atom(atom.name, element=atom.element, residue=res)
    allatom1 = md.Trajectory(allatom0.xyz, top, t.time, t.unitcell_lengths, t.unitcell_angles)
    allatom1.save_dcd(prot.path+'/allatom.dcd')
    allatom1[0].save_pdb(prot.path+'/allatom.pdb')
    print(prot_name,'has',allatom1.n_frames,'frames')

proteins = pd.read_pickle('proteins.pkl')

t0 = time.time()
reconstruct_pulchra(args.pulchra_path,args.name,proteins.loc[args.name],args.num_cpus)
print('Timing {:.3f}'.format(time.time()-t0))

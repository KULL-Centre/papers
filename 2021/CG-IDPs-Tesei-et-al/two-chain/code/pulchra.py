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
parser.add_argument('--ff',dest='ff',type=str)
parser.add_argument('--run',dest='run',type=int)
parser.add_argument('--pulchra',dest='pulchra_path',type=str)
args = parser.parse_args()

ray.init(num_cpus=args.num_cpus)

def fix_topology(dcd,pdb,n_chains):
    """
    Changes atom names to CA
    This function loops through all atoms/beads in 
    the topology changes the bead name to CA
    """
    # load trajectory
    t = md.load_dcd(dcd,pdb)
    # create new topology
    cgtop = md.Topology()  
    # add residues and CA atoms to new topology
    for _ in range(n_chains):
        cgchain = cgtop.add_chain()
        t_slice = t.atom_slice(range(0,int(t.n_atoms/n_chains)))
        for atom in t_slice.top.atoms:
            cgres = cgtop.add_residue(atom.name, cgchain, resSeq=atom.index+1)
            cgtop.add_atom('CA', element=md.element.carbon, residue=cgres)
    # create trajectory based on new topology
    traj = md.Trajectory(t.xyz, cgtop, t.time, t.unitcell_lengths, t.unitcell_angles)
    # fix one chain in box center and display closest periodic image of second chain
    traj = traj.image_molecules(inplace=False, anchor_molecules=[set(traj.top.chain(0).atoms)],
                                other_molecules=[set(traj.top.chain(1).atoms)], make_whole=False)
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

def reconstruct_pulchra(pulchra_path,prot_name,prot,num_cpus,n_chains,ff,run):
    """
    This function reconstructs an all-atom trajectory from a Calpha trajectory.
    For a trajectory with multiple chains, each chain will be reconstructed
    individually, and then the reconstructed trajectories will be joined.. 
    Input: trajectory .dcd and topology .pdb file.
    n_procs: number of processors to use in parallel.
    Return: A reconstructed mdtraj trajectory.
    """
    path = prot.path+ff+'/run{:d}'.format(run)
    # load trajectory data
    traj = fix_topology(path+'/{:s}.dcd'.format(prot_name),
                        path+'/{:s}.pdb'.format(prot_name),n_chains)
    traj = traj[66:]
    # loop to perform reconstruction on individual chains
    for q in range(n_chains):
        t = traj.atom_slice(range(q*(len(prot.fasta)),(len(prot.fasta))*(q+1)))
        name = path+'/0_chain{}.pdb'.format(q)
        t[0].save(name)
        subprocess.run([pulchra_path,name])    
        s = md.load_pdb(path+'/0_chain{}.rebuilt.pdb'.format(q))
        xyz = np.empty((0,s.n_atoms,3))
        xyz = np.append( xyz, s.xyz )
        num_cpus = num_cpus - 1
        # run pulchra on all frames
        for j in range(1, t.n_frames, num_cpus):
            n = j+num_cpus if j+num_cpus<t.n_frames else t.n_frames
            xyz = np.append( xyz, np.vstack(ray.get([run_pulchra.remote(path,pulchra_path,i,t[i]) for i in range(j,n)])) )
        # create trajectory with new coordinates
        allatom0 = md.Trajectory(xyz.reshape(t.n_frames,s.n_atoms,3), s.top, t.time, t.unitcell_lengths, t.unitcell_angles)
        # create topology with fixed residue indexes
        top = md.Topology()
        chain = top.add_chain()
        for residue in allatom0.top.residues:
            res = top.add_residue(residue.name, chain, resSeq=residue.index+1)
            for atom in residue.atoms:
                top.add_atom(atom.name, element=atom.element, residue=res)
        # use new topology to create the final all atom trajectory
        allatom1 = md.Trajectory(allatom0.xyz, top, t.time, t.unitcell_lengths, t.unitcell_angles)
        allatom1.save_dcd(path+'/allatom_chain{}.dcd'.format(q))
        allatom1[0].save_pdb(path+'/allatom_chain{}.pdb'.format(q))
        print(prot_name,ff,'has',allatom1.n_frames,'frames')
    # make stacked trajectory containing all-atom trajectories for both chains
    if n_chains == 2:
        t0 = md.load_dcd(path+'/allatom_chain0.dcd',path+'/allatom_chain0.pdb')
        t1 = md.load_dcd(path+'/allatom_chain1.dcd',path+'/allatom_chain1.pdb')
        stacked_traj = t0.stack(t1)
        stacked_traj.save_dcd(path+'/allatom.dcd')
        stacked_traj[0].save_pdb(path+'/allatom.pdb')
        # delete files with data on individual chains
        for chain_number in [0,1]:
            os.remove(path+'/0_chain{}.pdb'.format(chain_number))
            os.remove(path+'/0_chain{}.rebuilt.pdb'.format(chain_number))
            os.remove(path+'/allatom_chain{}.dcd'.format(chain_number))
            os.remove(path+'/allatom_chain{}.pdb'.format(chain_number))

proteins = pd.read_pickle('proteins.pkl')
n_chains = 2

t0 = time.time()
reconstruct_pulchra(args.pulchra_path,args.name,proteins.loc[args.name],args.num_cpus,n_chains,args.ff,args.run)
print('Timing {:.3f}'.format(time.time()-t0))

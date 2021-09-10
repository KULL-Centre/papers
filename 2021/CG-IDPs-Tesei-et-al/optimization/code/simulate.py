from analyse import *
import hoomd
import hoomd.md
import time
import os
import ray
import sys
import psutil
from argparse import ArgumentParser
from PeptideBuilder import Geometry
import PeptideBuilder
import Bio.PDB

parser = ArgumentParser()
parser.add_argument('--outdir',nargs='?',const='', type=str)
parser.add_argument('--num_cpus',dest='num_cpus',type=int)
args = parser.parse_args()

print(args.outdir,psutil.cpu_count())

replicas = lambda prot : int(np.ceil(len(prot.fasta)**2/3.8e4))
ray.init(num_cpus=args.num_cpus)

@ray.remote
def simulate(residues,name,prot,replica):
    Nsteps = np.ceil(3e7 / replicas(prot))
    hoomd.context.initialize("--mode=cpu --nthreads=1");
    hoomd.option.set_notice_level(0)
    hoomd.util.quiet_status()
    pairs, lj_eps, lj_lambda, lj_sigma, fasta, types, MWs = genParamsLJ(residues,name,prot)
    yukawa_eps, yukawa_kappa = genParamsDH(residues,name,prot)
    N = len(fasta)
    snapshot = hoomd.data.make_snapshot(N=N,
                                    box=hoomd.data.boxdim(Lx=200, Ly=200, Lz=200),
                                    particle_types=types,
                                    bond_types=['polymer']);

    geo = Geometry.geometry(prot.fasta[0])
    geo.phi = -120
    geo.psi_im1 = 150
    structure = PeptideBuilder.initialize_res(geo)
    for residue in prot.fasta[1:]:
        structure = PeptideBuilder.add_residue(structure, residue, geo.phi, geo.psi_im1)
    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    xyz = []
    for atom in out.structure.get_atoms():
        if atom.name == 'CA':
            xyz.append([atom.coord[0]/10.,atom.coord[1]/10.,atom.coord[2]/10.])
    snapshot.particles.position[:] = xyz;

    fst = ''.join(prot.fasta)
    ids = [types.index(a) for a in fasta]
    snapshot.particles.typeid[:] = ids
    snapshot.particles.mass[:] = [residues.loc[a].MW for a in prot.fasta]
    snapshot.particles.mass[0] += 2
    snapshot.particles.mass[-1] += 16
    snapshot.bonds.resize(N-1);
    snapshot.bonds.group[:] = [[i,i+1] for i in range(N-1)];
    snapshot.bonds.typeid[:] = [0] * (N-1)
    hoomd.init.read_snapshot(snapshot);

    kT = 8.3145*prot.temp*1e-3
    hb = hoomd.md.bond.harmonic();
    hb.bond_coeff.set('polymer', k=8033.0, r0=0.38);

    nl = hoomd.md.nlist.cell();

    lj1 = hoomd.md.pair.lj(r_cut=4.0, nlist=nl, name="1")
    lj2 = hoomd.md.pair.lj(r_cut=4.0, nlist=nl, name="2")
    yukawa = hoomd.md.pair.yukawa(r_cut=4.0, nlist=nl)
    for a,b in pairs:
        lj1.pair_coeff.set(a, b, epsilon=lj_eps*(1-lj_lambda.loc[a,b]), sigma=lj_sigma.loc[a,b], 
                    r_cut=np.power(2.,1./6)*lj_sigma.loc[a,b])
        lj2.pair_coeff.set(a, b, epsilon=lj_eps*lj_lambda.loc[a,b], sigma=lj_sigma.loc[a,b], r_cut=4.)
        yukawa.pair_coeff.set(a, b, epsilon=yukawa_eps.loc[a,b], kappa=yukawa_kappa, r_cut=4.);
    lj1.set_params(mode='shift')
    yukawa.set_params(mode='shift')
    nl.reset_exclusions(exclusions = ['bond'])

    integrator_mode = hoomd.md.integrate.mode_standard(dt=0.005);
    integrator = hoomd.md.integrate.langevin(group=hoomd.group.all(),kT=kT,seed=np.random.randint(100));
    for a,mw in zip(types,MWs):
        integrator.set_gamma(a, mw/100)
    nl.tune()
    hoomd.dump.gsd(prot.path+"/{:s}_{:d}.gsd".format(name,replica), period=2e3, group=hoomd.group.all(), overwrite=True);
    hoomd.run(Nsteps+2e6)

@ray.remote
def genDCD(residues,name,prot):
    """ 
    Generates coordinate and trajectory
    in convenient formats
    """
    top = md.Topology()
    chain = top.add_chain()
    for resname in prot.fasta:
        residue = top.add_residue(residues.loc[resname,'three'], chain)
        top.add_atom(residues.loc[resname,'three'], element=md.element.carbon, residue=residue)
    for i in range(len(prot.fasta)-1):
        top.add_bond(top.atom(i),top.atom(i+1))
    traj = md.load(prot.path+"/{:s}_0.gsd".format(name), top)[1000:]
    traj.top = top
    for replica in range(1,replicas(prot)):
        t = md.load(prot.path+"/{:s}_{:d}.gsd".format(name,replica), top)[1000:]
        t.top = top
        traj = md.join([traj,t])
        del t
    traj = traj.image_molecules(inplace=False, anchor_molecules=[set(traj.top.chain(0).atoms)], make_whole=True)
    traj.center_coordinates()
    traj.xyz += traj.unitcell_lengths[0,0]/2
    traj.save_dcd(prot.path+"/{:s}.dcd".format(name))
    traj[0].save_pdb(prot.path+"/{:s}.pdb".format(name))
    for replica in range(replicas(prot)):
        os.remove(prot.path+"/{:s}_{:d}.gsd".format(name,replica))

proteins = initProteins(args.outdir)
proteinsRgs = initProteinsRgs(args.outdir)

proteins.to_pickle('proteins.pkl')
proteinsRgs.to_pickle('proteinsRgs.pkl')
residues = pd.read_pickle('residues.pkl').set_index('one')

istraj = 0
for name, prot in proteins.iterrows():
    istraj += not os.path.isfile(prot.path+"/{:s}.dcd".format(name))
for name, prot in proteinsRgs.iterrows():
    istraj += not os.path.isfile(prot.path+"/{:s}.dcd".format(name))
print('Simulate?',bool(istraj))

allproteins = pd.concat((proteins,proteinsRgs),sort=True)

print('Number of simulations',len([r for name,prot in allproteins.iterrows() for r in range(replicas(prot))]))

if bool(istraj):
    t0 = time.time()
    ray.get([simulate.remote(residues,name,prot,r) for name,prot in allproteins.iterrows() for r in range(replicas(prot))])
    ray.get([genDCD.remote(residues,name,prot) for name,prot in allproteins.iterrows()])
    print('Timing Simulation {:.3f}'.format(time.time()-t0))

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
    pairs, lj_eps, lj_lambda, lj_sigma, fasta, types = genParamsLJ(residues,name,prot)
    yukawa_eps, yukawa_kappa = genParamsDH(residues,name,prot)
    N = len(fasta)
    snapshot = hoomd.data.make_snapshot(N=N,
                                    box=hoomd.data.boxdim(Lx=300, Ly=300, Lz=300),
                                    particle_types=types,
                                    bond_types=['polymer'],
                                    angle_types=['xgp', 'xpp', 'xgx', 'xpx', 'xxp', 'xxx'],
                                    dihedral_types=['gg','gp','pg','pp','gx','px','xg','xp','xx']);

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
    for i,idi in enumerate(ids):
        snapshot.particles.typeid[i] = idi
        snapshot.particles.mass[i] = residues.loc[types[idi]].MW
    snapshot.particles.mass[0] += 2
    snapshot.particles.mass[-1] += 16
    snapshot.bonds.resize(N-1);
    snapshot.bonds.group[:] = [[i,i+1] for i in range(N-1)];
    snapshot.bonds.typeid[:] = [0] * (N-1)
    snapshot.angles.resize(N-2);
    snapshot.angles.group[:] = [[i,i+1,i+2] for i in range(N-2)];
    snapshot.angles.typeid[:] = [0 if 'GP' == fst[i+1:i+3] else 1 if 'PP' == fst[i+1:i+3] else 2 if 'G' == fst[i+1] else 3 if 'P' == fst[i+1] else 4 if 'P' == fst[i+2] else 5 for i in range(N-2)]
    snapshot.dihedrals.resize(N-3);
    snapshot.dihedrals.group[:] = [[i,i+1,i+2,i+3] for i in range(N-3)];
    snapshot.dihedrals.typeid[:] = [0 if 'GG' == fst[i+1:i+3] else 1 if 'GP' == fst[i+1:i+3] else 2 if 'PG' == fst[i+1:i+3] else 3 if 'PP' == fst[i+1:i+3] else 4 if 'G' == fst[i+1] else 5 if 'P' == fst[i+1] else 6 if 'G' == fst[i+2] else 7 if 'P' == fst[i+2] else 8 for i in range(N-3)]
    hoomd.init.read_snapshot(snapshot);
    kT = 8.3145*prot.temp*1e-3
    hb = hoomd.md.bond.harmonic();
    hb.bond_coeff.set('polymer', k=8033.0, r0=0.38);

    atable = hoomd.md.angle.table(width=1000);
    dtable = hoomd.md.dihedral.table(width=2000);

    def angle_potential(x,p0,p1,p2,p3,p4,p5):
        V = p0*np.cos(x-p1) + p2*np.cos(2*x-p3) + p4*np.cos(3*x-p5)
        T = p0*np.sin(x-p1) + 2*p2*np.sin(2*x-p3) + 3*p4*np.sin(3*x-p5)
        return (V, T)

    def dihedral_potential(x,p0,p1,p2,p3):
        V = p0*np.cos(x-p1) + p2*np.cos(2*x-p3)
        T = p0*np.sin(x-p1) + 2*p2*np.sin(2*x-p3)
        return (V, T)

    paramsA = {'xgp': [-25544.630825164724, -4.250555383915038, 10662.64737109302, 4.067927233557355, 1912.3404546171032, 2.9701139887833987],
        'xpp': [-47614.27624849932, 1.9220723496497203, 19375.730077566997, 3.8544305004873745, 3344.1827255510507, 2.669005619847293],
        'xgx': [-63241.26002896861, -4.24306542628382, 26216.77249977118, 4.082304992696825, 4639.827915575768, 2.9877293946420234],
        'xpx': [141470.74112737243, -1.2041765233140187, 58055.23747945211, 3.8784923797633604, 10116.324086397397, 2.6859833895330176],
        'xxp': [-62523.829596244825, 1.9684990278030738, 26093.587759241862, 3.9373460546927133, 4680.183244125726, 2.7661636640387335],
        'xxx': [-52182.28992720884, 8.293194519031209, 21645.64446238669, 4.0210151766830275, 3836.5399042564536, 2.892670422941219]}

    paramsD = {'gg': [-0.05903219460599078, 2.034443935795617, -1.826128880993862, 3.031134993379038],
        'gp': [6.25179229981291, 1.76736372948421, 0.9229917713609358, 0.2640570072075076],
        'pg': [-0.9852586461250354, 2.3893576827589413, -2.4486200133954545, 1.5748394314489427],
        'pp': [7.751901593792672, 1.4932408389764515, -5.832925266105529, 2.048107053129765],
        'gx': [-1.340385093914546, 3.1120446057221, -2.09899975464143, 3.352246764263346],
        'px': [2.9313660569773634, 0.2619076846984134, -2.670778706669702, 1.4456750488732066],
        'xg': [-0.4834246580388573, 6.034895397876356, -2.067773258846332, 1.6618896773408387],
        'xp': [6.654937454248028, 0.8225740101666762, -2.6486010514757847, 1.4005268420293548],
        'xx': [1.6966621938382431, 0.02334206325518324, -1.4786982957993833, 1.5507097429670142]}

    for typeid,p in paramsA.items():
        atable.angle_coeff.set(typeid, func=angle_potential, coeff=dict(p0=p[0],p1=p[1],p2=p[2],p3=p[3],p4=p[4],p5=p[5]))

    for typeid,p in paramsD.items():
        dtable.dihedral_coeff.set(typeid, func=dihedral_potential, coeff=dict(p0=p[0],p1=p[1],p2=p[2],p3=p[3]))

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
    nl.reset_exclusions(exclusions = ['bond','angle'])

    integrator_mode = hoomd.md.integrate.mode_standard(dt=0.005);
    integrator = hoomd.md.integrate.langevin(group=hoomd.group.all(),kT=kT,seed=np.random.randint(100));
    for a in types:
        integrator.set_gamma(a, residues.loc[a].MW/100)
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

from analyse import *
import hoomd
import hoomd.md
import time
import os
import sys
from argparse import ArgumentParser
from mdtraj.utils.rotation import rotation_matrix_from_quaternion
from PeptideBuilder import Geometry
import PeptideBuilder
import Bio.PDB

parser = ArgumentParser()
parser.add_argument('--name',nargs='?',const='', type=str)
parser.add_argument('--temp',nargs='?',const='', type=int)
args = parser.parse_args()

print(hoomd.__file__)

def simulate(residues,name,prot,temp):
    residues = residues.set_index('one')
    hoomd.context.initialize("--mode=gpu");
    hoomd.option.set_notice_level(1)
    #hoomd.util.quiet_status()
    pairs, lj_eps, lj_lambda, lj_sigma, fasta, types, MWs = genParamsLJ(residues,name,prot)
    yukawa_eps, yukawa_kappa, charges = genParamsDH(residues,name,prot,temp)
    N = len(fasta)
    
    L = 15.
    if N > 400:
        L = 25.
        Lz = 300.
        Nsteps = 2e7
    elif N > 200:
        L = 17.
        Lz = 300.
        Nsteps = 6e7
    else:
        Lz = 10*L
        Nsteps = 6e7
    
    xy = np.empty(0)
    xy = np.append(xy,np.random.rand(2)*(L-2)-(L-2)/2).reshape((-1,2))
    for x,y in np.random.rand(1000,2)*(L-2)-(L-2)/2:
        x1 = x-L if x>0 else x+L
        y1 = y-L if y>0 else y+L
        if np.all(np.linalg.norm(xy-[x,y],axis=1)>.7):
            if np.all(np.linalg.norm(xy-[x1,y],axis=1)>.7):
                if np.all(np.linalg.norm(xy-[x,y1],axis=1)>.7):
                    xy = np.append(xy,[x,y]).reshape((-1,2))
        if xy.shape[0] == 100:
            break 
   
    n_chains = xy.shape[0]
    
    print('Number of chains',n_chains,N,'residues long')
           
    snapshot = hoomd.data.make_snapshot(N=N*n_chains,
                                box=hoomd.data.boxdim(Lx=L, Ly=L, Lz=Lz),
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
    xyz = np.array(xyz)
    v = xyz[-1] - xyz[0]
    u = np.array([0,0,1])
    a = np.cross(v,u) 
    a = a / np.linalg.norm(a,keepdims=True)
    b = np.arccos( np.dot(v,u) / np.linalg.norm(v) )
    quaternion = np.insert(np.sin(-b/2).reshape(-1,1)*a,0,np.cos(-b/2),axis=1)
    newxyz = xyz - np.mean(xyz,axis=0)
    newxyz = np.matmul(newxyz,rotation_matrix_from_quaternion(quaternion)) 
    xyz = np.array(newxyz[0])
    print(xyz[:,0].min(),xyz[:,0].max(),xy[:,0].min(),xy[:,0].max())
    print(xyz[:,1].min(),xyz[:,1].max(),xy[:,1].min(),xy[:,1].max())
    fst = ''.join(prot.fasta)

    snapshot.bonds.resize(n_chains*(N-1));

    for j,(x,y) in enumerate(xy):
        begin = j*N
        end = j*N+N
        
        snapshot.particles.position[begin:end] = [[xyz[i,0]+x,xyz[i,1]+y,xyz[i,2]] for i in range(N)];
        snapshot.particles.typeid[begin:end] = [types.index(a) for a in fasta]
        snapshot.particles.mass[begin:end] = [residues.loc[a].MW for a in prot.fasta]
        snapshot.particles.mass[begin] += 2
        snapshot.particles.mass[end-1] += 16
    
        snapshot.bonds.group[begin-j:end-j-1] = [[i,i+1] for i in range(begin,end-1)];
        snapshot.bonds.typeid[begin-j:end-j-1] = [0] * (N-1)

    hoomd.init.read_snapshot(snapshot);

    kT = 8.3145*temp*1e-3
    hb = hoomd.md.bond.harmonic();
    hb.bond_coeff.set('polymer', k=8033.0, r0=0.38);

    nl = hoomd.md.nlist.cell();

    lj1 = hoomd.md.pair.lj(r_cut=4.0, nlist=nl, name="one")
    lj2 = hoomd.md.pair.lj(r_cut=4.0, nlist=nl, name="two")
    yukawa = hoomd.md.pair.yukawa(r_cut=4.0, nlist=nl)
    for a,b in pairs:
        yukawa.pair_coeff.set(a, b, epsilon=yukawa_eps.loc[a,b], kappa=yukawa_kappa, r_cut=4.)
        lj1.pair_coeff.set(a, b, epsilon=lj_eps*(1-lj_lambda.loc[a,b]), sigma=lj_sigma.loc[a,b], 
                    r_cut=np.power(2.,1./6)*lj_sigma.loc[a,b])
        lj2.pair_coeff.set(a, b, epsilon=lj_eps*lj_lambda.loc[a,b], sigma=lj_sigma.loc[a,b], r_cut=4.)

    lj1.set_params(mode='shift')
    yukawa.set_params(mode='shift')
    nl.reset_exclusions(exclusions = ['bond'])

    walls = hoomd.md.wall.group()

    if N > 400:
        walls.add_plane((0,0,-50),(0,0,1))
        walls.add_plane((0,0,50),(0,0,-1))
    elif N > 200:
        walls.add_plane((0,0,-30),(0,0,1))
        walls.add_plane((0,0,30),(0,0,-1))
    else:
        walls.add_plane((0,0,-10),(0,0,1))
        walls.add_plane((0,0,10),(0,0,-1))
    gauss = hoomd.md.wall.gauss(walls,r_cut=4)
    gauss.force_coeff.set(types, sigma=1.0, epsilon=10.0)

    integrator_mode = hoomd.md.integrate.mode_standard(dt=0.005);   
    integrator = hoomd.md.integrate.langevin(group=hoomd.group.all(),kT=kT,seed=np.random.randint(100));

    for a,mw in zip(types,MWs):
        integrator.set_gamma(a, mw/100)
 
    hoomd.run(2e7)
    gauss.disable()

    hoomd.dump.gsd(filename=name+"/{:d}/{:s}.gsd".format(temp,name), period=5e4, group=hoomd.group.all(), overwrite=True);
    hoomd.dump.gsd(filename=name+"/{:d}/restart.gsd".format(temp), group=hoomd.group.all(), truncate=True, period=1e6, phase=0)

    hoomd.run(Nsteps)
    genDCD(residues,name,prot,temp,n_chains)

residues = pd.read_csv('residues.csv').set_index('three')
proteins = pd.read_pickle('proteins.pkl')
print(args.name,args.temp)
t0 = time.time()
simulate(residues,args.name,proteins.loc[args.name],args.temp)
print('Timing {:.3f}'.format(time.time()-t0))

from simtk import openmm, unit
from simtk.openmm import app
from simtk.openmm import XmlSerializer
from analyse import *
import time
import os
import sys
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--name',nargs='?',const='', type=str)
parser.add_argument('--temp',nargs='?',const='', type=int)
args = parser.parse_args()

def simulate(residues,name,prot,temp):
    residues = residues.set_index('one')

    lj_eps, fasta, types, MWs = genParamsLJ(residues,name,prot)
    yukawa_eps, yukawa_kappa = genParamsDH(residues,name,prot,temp)

    N = len(fasta)

    # set parameters
    L = 15.
    margin = 2
    if N > 400:
        L = 25.
        Lz = 300.
        margin = 8
        Nsteps = int(2e7)
    elif N > 200:
        L = 17.
        Lz = 300.
        margin = 4
        Nsteps = int(6e7)
    else:
        Lz = 10*L
        Nsteps = int(6e7)

    system = openmm.System()

    # set box vectors
    a = unit.Quantity(np.zeros([3]), unit.nanometers)
    a[0] = L * unit.nanometers
    b = unit.Quantity(np.zeros([3]), unit.nanometers)
    b[1] = L * unit.nanometers
    c = unit.Quantity(np.zeros([3]), unit.nanometers)
    c[2] = Lz * unit.nanometers
    system.setDefaultPeriodicBoxVectors(a, b, c)
    
    # initial config
    xy = np.empty(0)
    xy = np.append(xy,np.random.rand(2)*(L-margin)-(L-margin)/2).reshape((-1,2))
    for x,y in np.random.rand(1000,2)*(L-margin)-(L-margin)/2:
        x1 = x-L if x>0 else x+L
        y1 = y-L if y>0 else y+L
        if np.all(np.linalg.norm(xy-[x,y],axis=1)>.7):
            if np.all(np.linalg.norm(xy-[x1,y],axis=1)>.7):
                if np.all(np.linalg.norm(xy-[x,y1],axis=1)>.7):
                    xy = np.append(xy,[x,y]).reshape((-1,2))
        if xy.shape[0] == 100:
            break

    n_chains = xy.shape[0]

    top = md.Topology()
    pos = []
    for x,y in xy:
        chain = top.add_chain()
        pos.append([[x,y,Lz/2+(i-N/2.)*.38] for i in range(N)])
        for resname in fasta:
            residue = top.add_residue(resname, chain)
            top.add_atom(resname, element=md.element.carbon, residue=residue)
        for i in range(chain.n_atoms-1):
            top.add_bond(chain.atom(i),chain.atom(i+1))
    md.Trajectory(np.array(pos).reshape(n_chains*N,3), top, 0, [L,L,Lz], [90,90,90]).save_pdb(name+'/{:d}/top.pdb'.format(temp))

    pdb = app.pdbfile.PDBFile(name+'/{:d}/top.pdb'.format(temp))

    for _ in range(n_chains):
        system.addParticle((residues.loc[prot.fasta[0]].MW+2)*unit.amu)
        for a in prot.fasta[1:-1]:
            system.addParticle(residues.loc[a].MW*unit.amu) 
        system.addParticle((residues.loc[prot.fasta[-1]].MW+16)*unit.amu)

    hb = openmm.openmm.HarmonicBondForce()
    energy_expression = 'select(step(r-2^(1/6)*s),4*eps*l*((s/r)^12-(s/r)^6),4*eps*((s/r)^12-(s/r)^6)+eps*(1-l))'
    ah = openmm.openmm.CustomNonbondedForce(energy_expression+'; s=0.5*(s1+s2); l=0.5*(l1+l2)')
    yu = openmm.openmm.CustomNonbondedForce('q*(exp(-kappa*r)/r - exp(-kappa*4)/4); q=q1*q2')
    yu.addGlobalParameter('kappa',yukawa_kappa/unit.nanometer)
    yu.addPerParticleParameter('q')

    ah.addGlobalParameter('eps',lj_eps*unit.kilojoules_per_mole)
    ah.addPerParticleParameter('s')
    ah.addPerParticleParameter('l')
 
    for j in range(n_chains):
        begin = j*N
        end = j*N+N
       
        for a,e in zip(prot.fasta,yukawa_eps):
            yu.addParticle([e*unit.nanometer*unit.kilojoules_per_mole])
            ah.addParticle([residues.loc[a].sigmas*unit.nanometer, residues.loc[a].lambdas*unit.dimensionless])

        for i in range(begin,end-1):
            hb.addBond(i, i+1, 0.38*unit.nanometer, 8033.28*unit.kilojoules_per_mole/(unit.nanometer**2))
            yu.addExclusion(i, i+1)
            ah.addExclusion(i, i+1)

    yu.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
    ah.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
    hb.setUsesPeriodicBoundaryConditions(True)
    yu.setCutoffDistance(4*unit.nanometer)
    ah.setCutoffDistance(4*unit.nanometer)
 
    system.addForce(hb)
    system.addForce(yu)
    system.addForce(ah)

    serialized_system = XmlSerializer.serialize(system)
    outfile = open('system.xml','w')
    outfile.write(serialized_system)
    outfile.close()

    integrator = openmm.openmm.LangevinIntegrator(temp*unit.kelvin,0.01/unit.picosecond,0.005*unit.picosecond) # 322

    platform = openmm.Platform.getPlatformByName('CUDA')

    simulation = app.simulation.Simulation(pdb.topology, system, integrator, platform, dict(CudaPrecision='mixed'))

    check_point = name+'/{:d}/restart.chk'.format(temp)

    if os.path.isfile(check_point):
        simulation.loadCheckpoint(check_point)
        simulation.reporters.append(app.dcdreporter.DCDReporter(name+'/{:d}/{:s}.dcd'.format(temp,name),int(5e4),append=True))
    else:
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(app.dcdreporter.DCDReporter(name+'/{:d}/{:s}.dcd'.format(temp,name),int(5e4)))

    simulation.reporters.append(app.statedatareporter.StateDataReporter('{:s}_{:d}.log'.format(name,temp),100000,
             potentialEnergy=True,temperature=True,step=True,speed=True,elapsedTime=True,separator='\t'))

    simulation.runForClockTime(20*unit.hour, checkpointFile=check_point, checkpointInterval=5*unit.hour)

    simulation.saveCheckpoint(check_point)

    genDCD(residues,name,prot,temp,n_chains)

residues = pd.read_csv('residues.csv').set_index('three',drop=False)
proteins = pd.read_pickle('proteins.pkl')
print(args.name,args.temp)
t0 = time.time()
simulate(residues,args.name,proteins.loc[args.name],args.temp)
print('Timing {:.3f}'.format(time.time()-t0))

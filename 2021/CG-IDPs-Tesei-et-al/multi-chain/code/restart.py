from analyse import *
import hoomd
import hoomd.md
import time
import os
import sys
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--name',nargs='?',const='', type=str)
parser.add_argument('--temp',nargs='?',const='', type=int)
args = parser.parse_args()

print(hoomd.__file__)

def simulate(residues,name,prot,temp):
    residues = residues.set_index('one',drop=False)
    hoomd.context.initialize("--mode=gpu");
    hoomd.option.set_notice_level(1)
    hoomd.util.quiet_status()
    pairs, lj_eps, lj_lambda, lj_sigma, fasta, types, MWs = genParamsLJ(residues,name,prot)
    yukawa_eps, yukawa_kappa, charges = genParamsDH(residues,name,prot,temp)
    N = len(fasta)

    hoomd.init.read_gsd(name+"/{:d}/restart.gsd".format(temp,name));
    n_chains = int( len(hoomd.group.all()) / len(prot.fasta) )
 
    L = 15.
    if N > 200:
        L = 17.
        Lz = 300.
        Nsteps = 6e7
    else:
        Lz = 10*L
        Nsteps = 6e7

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

    integrator_mode = hoomd.md.integrate.mode_standard(dt=0.005);   
    integrator = hoomd.md.integrate.langevin(group=hoomd.group.all(),kT=kT,seed=np.random.randint(100));

    for a,mw in zip(types,MWs):
        integrator.set_gamma(a, mw/100)
   
    hoomd.dump.gsd(filename=name+"/{:d}/{:s}.gsd".format(temp,name), period=5e4, group=hoomd.group.all(), overwrite=False);
    hoomd.dump.gsd(filename=name+"/{:d}/restart.gsd".format(temp), group=hoomd.group.all(), truncate=True, period=1e6, phase=0)

    hoomd.run(Nsteps)
    genDCD(residues,name,prot,temp,n_chains)

residues = pd.read_csv('residues.csv').set_index('three',drop=False)
proteins = pd.read_pickle('proteins.pkl')
print(args.name,args.temp)
t0 = time.time()
simulate(residues,args.name,proteins.loc[args.name],args.temp)
print('Timing {:.3f}'.format(time.time()-t0))

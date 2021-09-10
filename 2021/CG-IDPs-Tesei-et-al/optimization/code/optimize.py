from analyse import *
import MDAnalysis
import time
import os
import glob
import sys
from DEERPREdict.PRE import PREpredict
from argparse import ArgumentParser
import psutil
import ray
import logging
import shutil
import h5py
import warnings
warnings.filterwarnings('ignore')

parser = ArgumentParser()
parser.add_argument('--log',dest='log_path',type=str,required=True)
parser.add_argument('--cycle',dest='cycle',type=str,required=True)
parser.add_argument('--num_cpus',dest='num_cpus',type=int)
args = parser.parse_args()

logging.basicConfig(filename=args.log_path+'/log',level=logging.INFO)

dp = 0.05
# set the confidence parameter
theta = 0.1

os.environ["NUMEXPR_MAX_THREADS"]="1"

proteins = pd.read_pickle('proteins.pkl')
proteinsRgs = pd.read_pickle('proteinsRgs.pkl')

for _, prot in proteins.iterrows():
    if not os.path.isdir(prot.path):
        os.mkdir(prot.path)
    if not os.path.isdir(prot.path+'/calcPREs'):
        os.mkdir(prot.path+'/calcPREs')

proc_PRE = [(label,name) for name,prot in proteins.iterrows() for label in prot.labels]
 
ray.init(num_cpus=args.num_cpus)

# function to calculate PRE data from all-atom trajectories using DEER-PREdict
@ray.remote
def evaluatePRE(n, label, name, prot):
    prefix = prot.path+'/calcPREs/res'
    filename = prefix+'-{:d}.pkl'.format(label)
    if isinstance(prot.weights, np.ndarray):
        u = MDAnalysis.Universe(prot.path+'/allatom.pdb')
        load_file = filename
    elif isinstance(prot.weights, bool):
        u = MDAnalysis.Universe(prot.path+'/allatom.pdb',prot.path+'/allatom.dcd')
        load_file = False
    else:
        raise ValueError('Weights argument is a '+str(type(prot.weights)))
    PRE = PREpredict(u, label, log_file = args.log_path+'/'+name+'/log', temperature = prot.temp, atom_selection = 'N', sigma_scaling = 1.0)
    PRE.run(output_prefix = prefix, weights = prot.weights, load_file = load_file, tau_t = 1e-10, tau_c = prot.tau_c*1e-09, r_2 = 10, wh = prot.wh)

# to speed-up reweighting, sum up the terms that are unaffected by changes in the lambda values and store them in a hdf5 data container
@ray.remote
def calcDistSums(df,name,prot):
    traj = md.load_dcd(prot.path+"/{:s}.dcd".format(name),prot.path+"/{:s}.pdb".format(name))
    
    _, _, _, lj_sigma, fasta, _, _ = genParamsLJ(df.set_index('one'),name,prot)
    pairs = traj.top.select_pairs('all','all')
    mask = np.abs(pairs[:,0]-pairs[:,1])>1 # exclude bonds
    pairs = pairs[mask]
    d = md.compute_distances(traj,pairs).astype(np.float32)
    d[d>4.] = np.inf # cutoff
    r = np.copy(d)
    n = np.zeros(r.shape,dtype=np.int8)
    pairs = np.array(list(itertools.combinations(fasta,2)))
    pairs = np.core.defchararray.add(pairs[:,0],pairs[:,1])
    pairs = pairs[mask] # exclude bonded
    for i,(a,b) in enumerate(pairs):
        mask = r[:,i]>np.power(2.,1./6)*lj_sigma.loc[a,b]
        r[:,i][mask] = np.inf # cutoff
        n[:,i][~mask] = 1
 
    unique_pairs = np.unique(pairs)

    # Store to hdf5 sum of d^-6 d^-12 for all unique aa pairs
    f = h5py.File(prot.path+'/dist_sums.hdf5', "w")
    f.create_dataset("unique", shape=unique_pairs.shape, data=unique_pairs.astype('S2'), compression="gzip")
    d12 = f.create_dataset("d12",
                (traj.n_frames, unique_pairs.size), fillvalue=0, compression="gzip")
    d6 = f.create_dataset("d6",
                (traj.n_frames, unique_pairs.size), fillvalue=0, compression="gzip")
    d12[:] = np.apply_along_axis(lambda x: pd.Series(index=pairs,data=x).groupby(level=0).sum(),
                          axis=1, arr=np.power(d,-12.))
    d6[:] = np.apply_along_axis(lambda x: pd.Series(index=pairs,data=x).groupby(level=0).sum(),
                          axis=1, arr=-np.power(d,-6.))
    r12 = f.create_dataset("r12",
                (traj.n_frames, unique_pairs.size), fillvalue=0, compression="gzip")
    r6 = f.create_dataset("r6",
                (traj.n_frames, unique_pairs.size), fillvalue=0, compression="gzip")
    ncut = f.create_dataset("ncut",
                (traj.n_frames, unique_pairs.size), fillvalue=0, compression="gzip")
    r12[:] = np.apply_along_axis(lambda x: pd.Series(index=pairs,data=x).groupby(level=0).sum(),
                          axis=1, arr=np.power(r,-12.))
    r6[:] = np.apply_along_axis(lambda x: pd.Series(index=pairs,data=x).groupby(level=0).sum(),
                          axis=1, arr=-np.power(r,-6.))
    ncut[:] = np.apply_along_axis(lambda x: pd.Series(index=pairs,data=x).groupby(level=0).sum(),
                          axis=1, arr=n)
    f.close()

# calculate per-frame AH energies (non-electrostatic, nonbonded) using the current lambda set
def calcLJenergy(df,name,prot):
    f = h5py.File(prot.path+'/dist_sums.hdf5', "r")
    unique = f.get('unique')[()].astype(str)

    _, lj_eps, lj_lambda, lj_sigma, _, _, _ = genParamsLJ(df.set_index('one'),name,prot)
    dflambda = lj_lambda.unstack()
    dflambda.index = dflambda.index.map('{0[0]}{0[1]}'.format)
    eps1 = lj_eps*(1-dflambda.loc[unique].values)
    eps2 = lj_eps*dflambda.loc[unique].values
    dfsigma = lj_sigma.unstack()
    dfsigma.index = dfsigma.index.map('{0[0]}{0[1]}'.format)
    sigma = dfsigma.loc[unique].values

    d6 = f.get('d6')
    d12 = f.get('d12')
    r6 = f.get('r6')
    r12 = f.get('r12')
    ncut = f.get('ncut')

    sigma6 = np.power(sigma,6)
    sigma12 = np.power(sigma6,2)
    lj_energy1 = (eps1*ncut+4*eps1*(sigma6*r6 + sigma12*r12)).sum(axis=1)
    lj_energy2 = (4*eps2*(sigma6*d6 + sigma12*d12)).sum(axis=1)
    f.close()
    return lj_energy1 + lj_energy2

# calculate weights based on the total per-frame AH energies calculated with the current lambda set and with the lambda set used for the simulations
@ray.remote
def calcWeights(df,name,prot):
    new_lj_energy = calcLJenergy(df,name,prot)

    lj_energy = np.loadtxt(prot.path+'/{:s}_LJenergy.dat'.format(name))
    kT = 8.3145*prot.temp*1e-3
    weights = np.exp((lj_energy-new_lj_energy)/kT)
    weights /= weights.sum()
    eff = np.exp(-np.sum(weights*np.log(weights*weights.size)))
    return (name,weights,eff)

# reweight the trajectory to estimate average observables for the current lambda set
def reweight(dp,df,dfprior,proteins,proteinsRgs):
    trial_proteins = proteins.copy()
    trial_proteinsRgs = proteinsRgs.copy()
    trial_df = df.copy()
    res_sel = np.random.choice(trial_df.index[:-2], 7, replace=False)
    trial_df.loc[res_sel,'lambdas'] += np.random.normal(0,dp,res_sel.size)
    f_out_of_01 = lambda df : df.loc[(df.lambdas<=0)|(df.lambdas>=1),'lambdas'].index
    f_out_of_prior = lambda df,dfprior : df.iloc[np.where(dfprior.lookup(df.one,df.lambdas)==0)[0]]['lambdas'].index

    out_of_01 = f_out_of_01(trial_df.iloc[:-2])
    trial_df.loc[out_of_01,'lambdas'] = df.loc[out_of_01,'lambdas']

    out_of_prior = f_out_of_prior(trial_df.iloc[:-2],dfprior)
    while (out_of_prior.size > 0):
        trial_df.loc[out_of_prior,'lambdas'] += np.random.normal(0,dp,out_of_prior.size)
        out_of_01 = f_out_of_01(trial_df.iloc[:-2])
        trial_df.loc[out_of_01,'lambdas'] = df.loc[out_of_01,'lambdas'] 
        out_of_prior = f_out_of_prior(trial_df.iloc[:-2],dfprior)

    prior = dfprior.lookup(trial_df.one[:-2],trial_df.lambdas[:-2])
    # calculate AH energies, weights and fraction of effective frames
    weights = ray.get([calcWeights.remote(trial_df,name,prot) for name,prot in pd.concat((trial_proteins,trial_proteinsRgs),sort=True).iterrows()])
    for name,w,eff in weights:
        if name in trial_proteins.index:
            trial_proteins.at[name,'weights'] = w
            trial_proteins.at[name,'eff'] = eff
        else:
            trial_proteinsRgs.at[name,'weights'] = w
            trial_proteinsRgs.at[name,'eff'] = eff
    # skip all sets of lambdas if the fraction of effective frames is too low
    if np.any(trial_proteins.eff < 0.3) or np.any(trial_proteinsRgs.eff < 0.3):
        return False, df, proteins, proteinsRgs, prior
    else:
        # calculate PREs and cost function
        ray.get([evaluatePRE.remote(n,label,name,trial_proteins.loc[name]) for n,(label,name) in enumerate(proc_PRE)])
        for name in trial_proteins.index:
            trial_proteins.at[name,'chi2_pre'] = calcChi2(trial_proteins.loc[name])
            rg, rh, rhkr, chi2_rh = calcRh(df,name,trial_proteins.loc[name])
            trial_proteins.at[name,'Rg'] = rg
            trial_proteins.at[name,'Rh'] = rh
            trial_proteins.at[name,'RhKR'] = rhkr
            trial_proteins.at[name,'chi2_rh'] = chi2_rh
        for name in trial_proteinsRgs.index:
            rg, chi2_rg = calcRg(df,name,trial_proteinsRgs.loc[name])
            trial_proteinsRgs.at[name,'Rg'] = rg
            trial_proteinsRgs.at[name,'chi2_rg'] = chi2_rg
        return True, trial_df, trial_proteins, trial_proteinsRgs, prior

df = pd.read_pickle('residues.pkl')

# set the initial lambda parameters to AVG
if args.cycle=='o1':
    df.lambdas = df.average

logging.info(df.lambdas-df.average)

# load the experimental PRE NMR data
for name in proteins.index:
    proteins.at[name,'expPREs'] = loadExpPREs(name,proteins.loc[name])

time0 = time.time()
ray.get([evaluatePRE.remote(n,label,name,proteins.loc[name]) for n,(label,name) in enumerate(proc_PRE)])
logging.info('Timing evaluatePRE {:.3f}'.format(time.time()-time0))

ray.get([calcDistSums.remote(df,name,prot) for name,prot in pd.concat((proteins,proteinsRgs),sort=True).iterrows()])

# initialize the protein DataFrames with the initial conformational properties
for name in proteins.index:
    lj_energy = calcLJenergy(df,name,proteins.loc[name])
    np.savetxt(proteins.loc[name].path+'/{:s}_LJenergy.dat'.format(name),lj_energy)
    tau_c, chi2_pre = optTauC(proteins.loc[name])
    proteins.at[name,'tau_c'] = tau_c
    proteins.at[name,'chi2_pre'] = chi2_pre
    rgarray, invrij, rg, rh, rhkr, chi2_rh = calcRh(df,name,proteins.loc[name])
    proteins.at[name,'rgarray'] = rgarray
    proteins.at[name,'invrij'] = invrij
    proteins.at[name,'Rg'] = rg
    proteins.at[name,'Rh'] = rh
    proteins.at[name,'RhKR'] = rhkr
    proteins.at[name,'chi2_rh'] = chi2_rh
    proteins.at[name,'initPREs'] = loadInitPREs(name,proteins.loc[name])
    if os.path.exists(proteins.loc[name].path+'/initPREs'):
        shutil.rmtree(proteins.loc[name].path+'/initPREs')
    shutil.copytree(proteins.loc[name].path+'/calcPREs',proteins.loc[name].path+'/initPREs')
proteins.to_pickle(args.cycle+'_init_proteins.pkl')

for name in proteinsRgs.index:
    lj_energy = calcLJenergy(df,name,proteinsRgs.loc[name])
    np.savetxt(proteinsRgs.loc[name].path+'/{:s}_LJenergy.dat'.format(name),lj_energy)
    rgarray, rg, chi2_rg = calcRg(df,name,proteinsRgs.loc[name])
    proteinsRgs.at[name,'rgarray'] = rgarray
    proteinsRgs.at[name,'Rg'] = rg
    proteinsRgs.at[name,'chi2_rg'] = chi2_rg
proteinsRgs.to_pickle(args.cycle+'_init_proteinsRgs.pkl')

logging.info('Initial Chi2 PRE {:.3f} +/- {:.3f}'.format(proteins.chi2_pre.mean(),proteins.chi2_pre.std()))
logging.info('Initial Chi2 Hydrodynamic Radius {:.3f} +/- {:.3f}'.format(proteins.chi2_rh.mean(),proteins.chi2_rh.std()))
logging.info('Initial Chi2 Gyration Radius {:.3f} +/- {:.3f}'.format(proteinsRgs.chi2_rg.mean(),proteinsRgs.chi2_rg.std()))

dfprior = pd.read_pickle('prior.pkl')
dfchi2 = pd.DataFrame(columns=['chi2_pre','chi2_rg','prior','w_aSyn','w_aSyn140','lambdas'])
dflambdas = pd.DataFrame(columns=['chi2_pre','chi2_rg','spearman','lambdas'])

# set the initial value for the control parameter of simulated annealing
if (args.cycle=='o1' or args.cycle=='o2'):
    kT0 = 2
else:
    kT0 = 0.1
logging.info('kT0 {:.1f}'.format(kT0))
kT = kT0

theta_prior = theta * np.log(dfprior.lookup(df.one[:-2],df.lambdas[:-2])).sum()

logging.info('Initial theta*prior {:.2f}'.format(theta_prior))

dfchi2.loc[0] = [proteins.chi2_pre.mean(),proteinsRgs.chi2_rg.mean(),theta_prior/theta,proteins.loc['aSyn','weights'],proteinsRgs.loc['aSyn140','weights'],df.lambdas[:-2]]

variants = ['A1','M12FP12Y','P7FM7Y','M9FP6Y','M8FP4Y','M9FP3Y','M10R','M6R','P2R','P7R','M3RP3K','M6RP6K','M10RP10K','M4D','P4D','P8D','P12D','P12E','P7KP12D','P7KP12Db','M12FP12YM10R','M10FP7RP12D']

# start simulated annealing
for k in range(2,5002):
    # interrupt simulated annealing when the control parameter drops below 1e-15
    if (kT<1e-15):
        logging.info('kT {:.2f}'.format(kT))
        break
    # scale the control parameter
    kT = kT * .99
    passed, trial_df, trial, trialRgs, prior = reweight(dp,df,dfprior,proteins,proteinsRgs)
    # passed is True if the effective fraction of frames exceeds 0.3
    if passed:
        theta_prior = theta * np.log(prior).sum()
        delta1 = trial.chi2_pre.mean() + trialRgs.chi2_rg.mean() - theta_prior
        delta2 = proteins.chi2_pre.mean() + proteinsRgs.chi2_rg.mean() - theta*dfchi2.iloc[-1]['prior']
        delta = delta1 - delta2
        # Metropolis criterion for accepting trial lambda parameters
        if ( np.exp(-delta/kT) > np.random.rand() ):
            proteins = trial.copy()
            proteinsRgs = trialRgs.copy()
            df = trial_df.copy()
            dfchi2.loc[k-1] = [trial.chi2_pre.mean(),trialRgs.chi2_rg.mean(),theta_prior/theta,trial.loc['aSyn','weights'],trialRgs.loc['aSyn140','weights'],df.lambdas[:-2]]
            if (trial.chi2_pre.mean() < 21) and (trialRgs.chi2_rg.mean() < 3):
                spearman = proteinsRgs.loc[variants,'expRg'].astype(float).corr(proteinsRgs.loc[variants,'Rg'].astype(float),method='spearman')
                dflambdas.loc[k-1] = [trial.chi2_pre,trialRgs.chi2_rg,spearman,df.lambdas]
            logging.info('Iter {:d}, PRE {:.2f}, Chi2 Rh {:.2f}, Chi2 Rg {:.2f}, Chi2 PRE Std {:.2f}, Chi2 Rg Std {:.2f}, theta*prior {:.2f}'.format(k-1,trial.chi2_pre.mean(),trial.chi2_rh.mean(),trialRgs.chi2_rg.mean(),trial.chi2_pre.std(),trialRgs.chi2_rg.std(),theta_prior))

dfchi2.to_pickle(args.cycle+'_chi2.pkl')
dflambdas.to_pickle(args.cycle+'_lambdas.pkl')
proteins.to_pickle('proteins.pkl')
proteins.to_pickle(args.cycle+'_proteins.pkl')
proteinsRgs.to_pickle('proteinsRgs.pkl')
proteinsRgs.to_pickle(args.cycle+'_proteinsRgs.pkl')
df.to_pickle('residues.pkl')
df.to_pickle(args.cycle+'_residues.pkl')

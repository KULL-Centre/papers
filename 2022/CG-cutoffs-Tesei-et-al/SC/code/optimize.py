from analyse import *
import MDAnalysis
import time
import os
import glob
import sys
from DEERPREdict.PRE import PREpredict
from argparse import ArgumentParser
from sklearn.neighbors import KernelDensity
import ray
import logging
import shutil
import warnings

warnings.filterwarnings('ignore')

parser = ArgumentParser()
parser.add_argument('--log',dest='log_path',type=str,required=True)
parser.add_argument('--cycle',dest='cycle',type=int,required=True)
parser.add_argument('--num_cpus',dest='num_cpus',type=int)
parser.add_argument('--cutoff',dest='cutoff',type=float)
args = parser.parse_args()

logging.basicConfig(filename=args.log_path+'/log',level=logging.INFO)

dp = 0.05
theta = .05
eta = .1
xi_0 = .1
rc = args.cutoff

os.environ["NUMEXPR_MAX_THREADS"]="1"

proteins = initProteins(args.cycle)
proteinsRgs = initProteinsRgs(args.cycle)

proteins.to_pickle('proteins.pkl')
proteinsRgs.to_pickle('proteinsRgs.pkl')

proteins = pd.read_pickle('proteins.pkl').astype(object)
proteinsRgs = pd.read_pickle('proteinsRgs.pkl').astype(object)

for _, prot in proteins.iterrows():
    if not os.path.isdir(prot.path):
        os.mkdir(prot.path)
    if not os.path.isdir(prot.path+'/calcPREs'):
        os.mkdir(prot.path+'/calcPREs')

proc_PRE = [(label,name) for name,prot in proteins.iterrows() for label in prot.labels]

ray.init(num_cpus=args.num_cpus)

@ray.remote(num_cpus=1)
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

@ray.remote(num_cpus=1)
def calcDistSums(df,name,prot):
    if not os.path.isfile(prot.path+'/energy_sums_2.npy'):
        traj = md.load_dcd(prot.path+"/{:s}.dcd".format(name),prot.path+"/{:s}.pdb".format(name))
        fasta = [res.name for res in traj.top.atoms]
        pairs = traj.top.select_pairs('all','all')
        mask = np.abs(pairs[:,0]-pairs[:,1])>1 # exclude bonds
        pairs = pairs[mask]
        d = md.compute_distances(traj,pairs)
        d[d>rc] = np.inf # cutoff
        r = np.copy(d)
        n1 = np.zeros(r.shape,dtype=np.int8)
        n2 = np.zeros(r.shape,dtype=np.int8)
        pairs = np.array(list(itertools.combinations(fasta,2)))
        pairs = pairs[mask] # exclude bonded
        sigmas = 0.5*(df.loc[pairs[:,0]].sigmas.values+df.loc[pairs[:,1]].sigmas.values)
        for i,sigma in enumerate(sigmas):
            mask = r[:,i]>np.power(2.,1./6)*sigma
            r[:,i][mask] = np.inf # cutoff
            n1[:,i][~mask] = 1
            n2[:,i][np.isfinite(d[:,i])] = 1

        unique_pairs = np.unique(pairs,axis=0)
        pairs = np.core.defchararray.add(pairs[:,0],pairs[:,1])

        d12 = np.apply_along_axis(lambda x: pd.Series(index=pairs,data=x).groupby(level=0).sum(),
                              axis=1, arr=np.power(d,-12.))
        d6 = np.apply_along_axis(lambda x: pd.Series(index=pairs,data=x).groupby(level=0).sum(),
                              axis=1, arr=-np.power(d,-6.))
        r12 = np.apply_along_axis(lambda x: pd.Series(index=pairs,data=x).groupby(level=0).sum(),axis=1, arr=np.power(r,-12.))
        r6 = np.apply_along_axis(lambda x: pd.Series(index=pairs,data=x).groupby(level=0).sum(),axis=1, arr=-np.power(r,-6.))
        ncut1 = np.apply_along_axis(lambda x: pd.Series(index=pairs,data=x).groupby(level=0).sum(),axis=1, arr=n1)
        ncut2 = np.apply_along_axis(lambda x: pd.Series(index=pairs,data=x).groupby(level=0).sum(),axis=1, arr=n2)
        sigmas = 0.5*(df.loc[unique_pairs[:,0]].sigmas.values+df.loc[unique_pairs[:,1]].sigmas.values)
        sigmas6 = np.power(sigmas,6)
        sigmas12 = np.power(sigmas6,2)
        eps = 0.2*4.184
        term_1 = eps*(ncut1+4*(sigmas6*r6 + sigmas12*r12))
        term_2 = sigmas6*(d6+ncut2/np.power(rc,6)) + sigmas12*(d12+ncut2/np.power(rc,12))
        term_2 = 4*eps*term_2 - term_1
        np.save(prot.path+'/energy_sums_1.npy',term_1.sum(axis=1))
        np.save(prot.path+'/energy_sums_2.npy',term_2)
        np.save(prot.path+'/unique_pairs.npy',unique_pairs)

def calcAHenergy(df,name,prot):
    term_1 = np.load(prot.path+'/energy_sums_1.npy')
    term_2 = np.load(prot.path+'/energy_sums_2.npy')
    unique_pairs = np.load(prot.path+'/unique_pairs.npy')
    lambdas = 0.5*(df.loc[unique_pairs[:,0]].lambdas.values+df.loc[unique_pairs[:,1]].lambdas.values)
    return term_1+np.nansum(lambdas*term_2,axis=1)

@ray.remote(num_cpus=1)
def calcWeights(df,name,prot):
    new_ah_energy = calcAHenergy(df,name,prot)
    ah_energy = np.load(prot.path+'/{:s}_AHenergy.npy'.format(name))
    kT = 8.3145*prot.temp*1e-3
    weights = np.exp((ah_energy-new_ah_energy)/kT)
    weights /= weights.sum()
    eff = np.exp(-np.sum(weights*np.log(weights*weights.size)))
    return name,weights,eff

def reweight(dp,df,proteins,proteinsRgs):
    trial_proteins = proteins.copy()
    trial_proteinsRgs = proteinsRgs.copy()
    trial_df = df.copy()
    res_sel = np.random.choice(trial_df.index, 5, replace=False)
    trial_df.loc[res_sel,'lambdas'] += np.random.normal(0,dp,res_sel.size)
    f_out_of_01 = lambda df : df.loc[(df.lambdas<=0)|(df.lambdas>1),'lambdas'].index
    out_of_01 = f_out_of_01(trial_df)
    trial_df.loc[out_of_01,'lambdas'] = df.loc[out_of_01,'lambdas']

    # calculate AH energies, weights and fraction of effective frames
    weights = ray.get([calcWeights.remote(trial_df,name,prot) for name,prot in pd.concat((trial_proteins,trial_proteinsRgs),sort=True).iterrows()])
    for name,w,eff in weights:
        if eff < 0.6:
            return False, df, proteins, proteinsRgs
        if name in trial_proteins.index:
            trial_proteins.at[name,'weights'] = w
            trial_proteins.at[name,'eff'] = eff
        else:
            trial_proteinsRgs.at[name,'weights'] = w
            trial_proteinsRgs.at[name,'eff'] = eff

    # calculate PREs and cost function
    ray.get([evaluatePRE.remote(n,label,name,trial_proteins.loc[name]) for n,(label,name) in enumerate(proc_PRE)])
    for name in trial_proteins.index:
        trial_proteins.at[name,'chi2_pre'] = calcChi2(trial_proteins.loc[name])
    for name in trial_proteinsRgs.index:
        rg, chi2_rg = reweightRg(df,name,trial_proteinsRgs.loc[name])
        trial_proteinsRgs.at[name,'Rg'] = rg
        trial_proteinsRgs.at[name,'chi2_rg'] = chi2_rg
    return True, trial_df, trial_proteins, trial_proteinsRgs

df = pd.read_csv('residues.csv').set_index('three')
logging.info(df.lambdas-df.M1)

for name in proteins.index:
    proteins.at[name,'expPREs'] = loadExpPREs(name,proteins.loc[name])

time0 = time.time()
ray.get([evaluatePRE.remote(n,label,name,proteins.loc[name]) for n,(label,name) in enumerate(proc_PRE)])
logging.info('Timing evaluatePRE {:.3f}'.format(time.time()-time0))

time0 = time.time()
ray.get([calcDistSums.remote(df,name,prot) for name,prot in pd.concat((proteins,proteinsRgs),sort=True).iloc[:10].iterrows()])
ray.get([calcDistSums.remote(df,name,prot) for name,prot in pd.concat((proteins,proteinsRgs),sort=True).iloc[10:20].iterrows()])
ray.get([calcDistSums.remote(df,name,prot) for name,prot in pd.concat((proteins,proteinsRgs),sort=True).iloc[20:30].iterrows()])
ray.get([calcDistSums.remote(df,name,prot) for name,prot in pd.concat((proteins,proteinsRgs),sort=True).iloc[30:40].iterrows()])
ray.get([calcDistSums.remote(df,name,prot) for name,prot in pd.concat((proteins,proteinsRgs),sort=True).iloc[40:50].iterrows()])
ray.get([calcDistSums.remote(df,name,prot) for name,prot in pd.concat((proteins,proteinsRgs),sort=True).iloc[50:].iterrows()])
logging.info('Timing calcDistSums {:.3f}'.format(time.time()-time0))

for name in proteins.index:
    np.save(proteins.loc[name].path+'/{:s}_AHenergy.npy'.format(name),
        calcAHenergy(df,name,proteins.loc[name]))
    tau_c, chi2_pre = optTauC(proteins.loc[name])
    proteins.at[name,'tau_c'] = tau_c
    proteins.at[name,'chi2_pre'] = chi2_pre
    proteins.at[name,'initPREs'] = loadInitPREs(name,proteins.loc[name])
    if os.path.exists(proteins.loc[name].path+'/initPREs'):
        shutil.rmtree(proteins.loc[name].path+'/initPREs')
    shutil.copytree(proteins.loc[name].path+'/calcPREs',proteins.loc[name].path+'/initPREs')
proteins.to_pickle(str(args.cycle)+'_init_proteins.pkl')

for name in proteinsRgs.index:
    np.save(proteinsRgs.loc[name].path+'/{:s}_AHenergy.npy'.format(name),
        calcAHenergy(df,name,proteinsRgs.loc[name]))
    rgarray, rg, chi2_rg = calcRg(df,name,proteinsRgs.loc[name])
    proteinsRgs.at[name,'rgarray'] = rgarray
    proteinsRgs.at[name,'Rg'] = rg
    proteinsRgs.at[name,'chi2_rg'] = chi2_rg
proteinsRgs.to_pickle(str(args.cycle)+'_init_proteinsRgs.pkl')

logging.info('Initial Chi2 PRE {:.3f} +/- {:.3f}'.format(proteins.chi2_pre.mean(),proteins.chi2_pre.std()))
logging.info('Initial Chi2 Gyration Radius {:.3f} +/- {:.3f}'.format(proteinsRgs.chi2_rg.mean(),proteinsRgs.chi2_rg.std()))

selHPS = pd.read_csv('selHPS.csv',index_col=0)

kde = KernelDensity(kernel='gaussian',bandwidth=0.05).fit(selHPS.T.values)
theta_prior = theta * kde.score_samples(df.lambdas.values.reshape(1, -1))[0]

xi = xi_0

variants = ['A1','M12FP12Y','P7FM7Y','M9FP6Y','M8FP4Y','M9FP3Y','M10R','M6R','P2R','P7R','M3RP3K','M6RP6K','M10RP10K','M4D','P4D','P8D','P12D','P12E','P7KP12D','P7KP12Db','M12FP12YM10R','M10FP7RP12D']

logging.info('Initial theta*prior {:.2f}'.format(theta_prior))
logging.info('theta {:.2f}'.format(theta))
logging.info('xi {:g}'.format(xi))

dfchi2 = pd.DataFrame(columns=['chi2_pre','chi2_rg','theta_prior','lambdas','xi'])
dfchi2.loc[0] = [proteins.chi2_pre.mean(),proteinsRgs.chi2_rg.mean(),theta_prior,df.lambdas,xi]

time0 = time.time()

micro_cycle = 0

logging.info(df.lambdas)
df['lambdas_{:d}'.format(args.cycle-1)] = df.lambdas

for k in range(2,200000):

    if (xi<1e-8):
        xi = xi_0
        micro_cycle += 1
        if (micro_cycle==10):
            logging.info('xi {:g}'.format(xi))
            break

    xi = xi * .99
    passed, trial_df, trial, trialRgs = reweight(dp,df,proteins,proteinsRgs)
    if passed:
        theta_prior = theta * kde.score_samples(trial_df.lambdas.values.reshape(1, -1))[0]
        delta1 = eta*trial.chi2_pre.mean() + trialRgs.chi2_rg.mean() - theta_prior
        delta2 = eta*proteins.chi2_pre.mean() + proteinsRgs.chi2_rg.mean() - dfchi2.iloc[-1]['theta_prior']
        delta = delta1 - delta2
        if ( np.exp(-delta/xi) > np.random.rand() ):
            proteins = trial.copy()
            proteinsRgs = trialRgs.copy()
            df = trial_df.copy()
            dfchi2.loc[k-1] = [trial.chi2_pre.mean(),trialRgs.chi2_rg.mean(),theta_prior,df.lambdas,xi]
            logging.info('Acc Iter {:d}, micro cycle {:d}, xi {:g}, Chi2 PRE {:.2f}, Chi2 Rg {:.2f}, theta*prior {:.2f}'.format(k-1,micro_cycle,xi,trial.chi2_pre.mean(),trialRgs.chi2_rg.mean(),theta_prior))

logging.info('Timing Reweighting {:.3f}'.format(time.time()-time0))
logging.info('Theta {:.3f}'.format(theta))

dfchi2['cost'] = dfchi2.chi2_rg + eta*dfchi2.chi2_pre - dfchi2.theta_prior
dfchi2.to_pickle(str(args.cycle)+'_chi2.pkl')
df.lambdas = dfchi2.loc[dfchi2.cost.idxmin()].lambdas
logging.info(df.lambdas)
df['lambdas_{:d}'.format(args.cycle)] = df.lambdas
proteins.to_pickle('proteins.pkl')
proteinsRgs.to_pickle('proteinsRgs.pkl')
df.to_csv('residues.csv')
logging.info('Cost at 0: {:.2f}'.format(dfchi2.loc[0].cost))
logging.info('Min Cost at {:d}: {:.2f}'.format(dfchi2.cost.idxmin(),dfchi2.cost.min()))


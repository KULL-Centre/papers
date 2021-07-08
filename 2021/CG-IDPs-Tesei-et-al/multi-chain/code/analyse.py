import pandas as pd
import numpy as np
import mdtraj as md
import itertools
import os
import MDAnalysis
from MDAnalysis import transformations

def initProteins():
    proteins = pd.DataFrame(index=['ht40','FUS','M10R','P7FM7Y','M8FP4Y','M4D','P4D','M6R','P8D','M9FP6Y','M12FP12Y','M9FP3Y','A1','AroP','AroM','AroMM','Ddx4WT','Ddx4CS','Ddx4RK','Ddx4FA','A2','A2NS','M10FP7RP12D','M6RP6K','P7R','M3RP3K','P7KP12D','P2R','P12D','A1NLS'],
            columns=['eps_factor','temp','expRg','expRgErr','Rg','rgarray','eff','chi2_rg','weights','pH','ionic','fasta'])
    fasta_ht40 = """MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEP
GSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKK
AKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREP
KKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGS
VQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTD
HGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL""".replace('\n', '') 
    fasta_FUS = """MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYG
QSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQS
SSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNS""".replace('\n', '')
    fasta_Ddx4WT = """MGDEDWEAEINPHMSSYVPIFEKDRYSGENGDNFNRTPASSSEMDDGPSR
RDHFMKSGFASGRNFGNRDAGECNKRDNTSTMGGFGVGKSFGNRGFSNSR
FEDGDSSGFWRESSNDCEDNPTRNRGFSKRGGYRDGNNSEASGPYRRGGR
GSFRGCRGGFGLGSPNNDLDPDECMQRTGGLFGSRRPVLSGTGNGDTSQS
RSGSGSERGGYKGLNEEVITGSGKNSWKSEAEGGES""".replace('\n', '')
    fasta_Ddx4CS = """MGDRDWRAEINPHMSSYVPIFEKDRYSGENGRNFNDTPASSSEMRDGPSE
RDHFMKSGFASGDNFGNRDAGKCNERDNTSTMGGFGVGKSFGNEGFSNSR
FERGDSSGFWRESSNDCRDNPTRNDGFSDRGGYEKGNNSEASGPYERGGR
GSFDGCRGGFGLGSPNNRLDPRECMQRTGGLFGSDRPVLSGTGNGDTSQS
RSGSGSERGGYKGLNEKVITGSGENSWKSEARGGES""".replace('\n', '')
    fasta_Ddx4RK = """MGDEDWEAEINPHMSSYVPIFEKDKYSGENGDNFNKTPASSSEMDDGPSK
KDHFMKSGFASGKNFGNKDAGECNKKDNTSTMGGFGVGKSFGNKGFSNSK
FEDGDSSGFWKESSNDCEDNPTKNKGFSKKGGYKDGNNSEASGPYKKGGK
GSFKGCKGGFGLGSPNNDLDPDECMQKTGGLFGSKKPVLSGTGNGDTSQS
KSGSGSEKGGYKGLNEEVITGSGKNSWKSEAEGGES""".replace('\n', '')
    fasta_Ddx4FA = """MGDEDWEAEINPHMSSYVPIAEKDRYSGENGDNANRTPASSSEMDDGPSR
RDHAMKSGAASGRNAGNRDAGECNKRDNTSTMGGAGVGKSAGNRGASNSR
AEDGDSSGAWRESSNDCEDNPTRNRGASKRGGYRDGNNSEASGPYRRGGR
GSARGCRGGAGLGSPNNDLDPDECMQRTGGLAGSRRPVLSGTGNGDTSQS
RSGSGSERGGYKGLNEEVITGSGKNSWKSEAEGGES""".replace('\n', '')
    fasta_A2NS = """GHMGRGGSFGFGDSRGGGGSFGPGPGSSFRGGSDGYGSGR
GFGDGYSGYGGGPGGGSFGGSPGYGGGRGGYGGGGPGYGSQGGGYGGGYDSYGGGSYGSG
SYSDFGSYSQQPSSYGPMKSGSFGGSRSMGGPYGGGSYGPGGSGGSGGYGGRSRY""".replace('\n', '')
    fasta_M10FP7RP12D = """GSMASADSSQRDRDDRGNFGDGRGGGGGGNDNFGRGGNGSDRGGGGGSRGDGRYGGDGDRYNGGGNDGRNGGGGGSYNDGGNYNNQ
SSNGDPMKGGNGRDRSSGPYDRGGQYGAKPRNQGGYGGSSSSRSYGSDRRG""".replace('\n', '')
    fasta_P12D = """GSMASADSSQRDRDDSGNFGDGRGGGFGGNDNFGRGGNFSDRGGFGGSRGDGGYGGDGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFDPMKGGNFGDRSSGPYDGGGQYFAKPRNQGGYGGSSSSSSYGSDRRF""".replace('\n', '')
    fasta_M6RP6K = """GSMASASSSQKGKSGSGNFGGGRGGGFGGNDNFGKGGNFSGRGGFGGSKGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFGGKSSGGSGGGGQYFAKPRNQGGYGGSSSSSSYGSGRKF""".replace('\n', '')
    fasta_P7R = """GSMASASSSQRGRSGRGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGRYGGSGDRYNGFGNDGRNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFRGRSSGPYGRGGQYFAKPRNQGGYGGSSSSRSYGSGRRF""".replace('\n', '')
    fasta_M4D = """GSMASASSSQRGRSGSGNFGGGRGGGFGGNGNFGRGGNFSGRGGFGGSRGGGGYGGSGGGYNGFGNSGSNFGGGGSYNGFGNYNNQ
SSNFGPMKGGNFGGRSSGPYGGGGQYFAKPRNQGGYGGSSSSSSYGSGRRF""".replace('\n', '')
    fasta_P4D = """GSMASASSSQRDRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGDFGGSRGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFGGRSSDPYGGGGQYFAKPRNQGGYGGSSSSSSYDSGRRF""".replace('\n', '')
    fasta_A1 = """GSMASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFGGRSSGGSGGGGQYFAKPRNQGGYGGSSSSSSYGSGRRF""".replace('\n', '')
    fasta_A1NLS = """GSMASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFGGRSSGPYGGGGQYFAKPRNQGGYGGSSSSSSYGSGRRF""".replace('\n', '')
    fasta_P8D = """GSMASASSSQRDRSGSGNFGGGRDGGFGGNDNFGRGDNFSGRGDFGGSRDGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFGGRSSDPYGGGGQYFAKPRNQDGYGGSSSSSSYDSGRRF""".replace('\n', '')
    fasta_M6R = """GSMASASSSQGGRSGSGNFGGGRGGGFGGNDNFGGGGNFSGSGGFGGSRGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFGGSSSGPYGGGGQYFAKPGNQGGYGGSSSSSSYGSGGRF""".replace('\n', '')
    fasta_M9FP6Y = """GSMASASSSQRGRSGSGNFGGGRGGGYGGNDNYGRGGNYSGRGGFGGSRGGGGYGGSGDGYNGGGNDGSNYGGGGSYNDSGNYNNQ
SSNFGPMKGGNYGGRSSGGSGGGGQYGAKPRNQGGYGGSSSSSSYGSGRRY""".replace('\n', '')
    fasta_M8FP4Y = """GSMASASSSQRGRSGSGNFGGGRGGGYGGNDNGGRGGNYSGRGGFGGSRGGGGYGGSGDGYNGGGNDGSNYGGGGSYNDSGNYNNQ
SSNFGPMKGGNYGGRSSGGSGGGGQYGAKPRNQGGYGGSSSSSSYGSGRRF""".replace('\n', '')
    fasta_P7FM7Y = """GSMASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGFGGSGDGFNGFGNDGSNFGGGGSFNDFGNFNNQ
SSNFGPMKGGNFGGRSSGGSGGGGQFFAKPRNQGGFGGSSSSSSFGSGRRF""".replace('\n', '')
    fasta_M10R = """GSMASASSSQGGSSGSGNFGGGGGGGFGGNDNFGGGGNFSGSGGFGGSGGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFGGSSSGPYGGGGQYFAKPGNQGGYGGSSSSSSYGSGGGF""".replace('\n', '')
    fasta_P2R = """GSMASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGYGGSGDGYNGFRNDGSNFGGGGRYNDFGNYNNQ
SSNFGPMKGGNFGGRSSGPYGGGGQYFAKPRNQGGYGGSSSSSSYGSGRRF""".replace('\n', '')
    fasta_P7KP12D = """GSMASADSSQRDRDDKGNFGDGRGGGFGGNDNFGRGGNFSDRGGFGGSRGDGKYGGDGDKYNGFGNDGKNFGGGGSYNDFGNYNNQ
SSNFDPMKGGNFKDRSSGPYDKGGQYFAKPRNQGGYGGSSSSKSYGSDRRF""".replace('\n', '')
    fasta_M3RP3K = """GSMASASSSQRGKSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSKGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFGGRSSGGSGGGGQYFAKPRNQGGYGGSSSSSSYGSGRKF""".replace('\n', '')
    fasta_M12FP12Y = """GSMASASSSQRGRSGSGNYGGGRGGGYGGNDNYGRGGNYSGRGGYGGSRGGGGYGGSGDGYNGYGNDGSNYGGGGSYNDYGNYNNQ
SSNYGPMKGGNYGGRSSGGSGGGGQYYAKPRNQGGYGGSSSSSSYGSGRRY""".replace('\n', '')
    fasta_M9FP3Y = """GSMASASSSQRGRSGSGNFGGGRGGGYGGNDNGGRGGNYSGRGGFGGSRGGGGYGGSGDGYNGGGNDGSNYGGGGSYNDSGNGNNQ
SSNFGPMKGGNYGGRSSGGSGGGGQYGAKPRNQGGYGGSSSSSSYGSGRRS""".replace('\n', '')
        
    proteins.loc['FUS'] = dict(eps_factor=0.2,temp=277,expRg=3.32,expRgErr=0.07,pH=5.5,fasta=list(fasta_FUS),ionic=0.15)
    proteins.loc['A1'] = dict(eps_factor=0.2,temp=296,expRg=2.61,expRgErr=0.11,pH=7.0,fasta=list(fasta_A1),ionic=0.15)
    proteins.loc['A1NLS'] = dict(eps_factor=0.2,temp=296,expRg=2.61,expRgErr=0.11,pH=7.0,fasta=list(fasta_A1NLS),ionic=0.15)
    proteins.loc['P12D'] = dict(eps_factor=0.2,temp=296,expRg=2.801,expRgErr=0.11,pH=7.0,fasta=list(fasta_P12D),ionic=0.15)
    proteins.loc['M10FP7RP12D'] = dict(eps_factor=0.2,temp=296,expRg=2.86,expRgErr=0.11,pH=7.0,fasta=list(fasta_M10FP7RP12D),ionic=0.15)
    proteins.loc['M6RP6K'] = dict(eps_factor=0.2,temp=296,expRg=2.787,expRgErr=0.11,pH=7.0,fasta=list(fasta_M6RP6K),ionic=0.15)
    proteins.loc['M4D'] = dict(eps_factor=0.2,temp=296,expRg=2.787,expRgErr=0.11,pH=7.0,fasta=list(fasta_M4D),ionic=0.15)
    proteins.loc['P4D'] = dict(eps_factor=0.2,temp=296,expRg=2.787,expRgErr=0.11,pH=7.0,fasta=list(fasta_P4D),ionic=0.15)
    proteins.loc['M6R'] = dict(eps_factor=0.2,temp=296,expRg=2.787,expRgErr=0.11,pH=7.0,fasta=list(fasta_M6R),ionic=0.15)
    proteins.loc['M10R'] = dict(eps_factor=0.2,temp=296,expRg=2.787,expRgErr=0.11,pH=7.0,fasta=list(fasta_M10R),ionic=0.15)
    proteins.loc['P8D'] = dict(eps_factor=0.2,temp=296,expRg=2.787,expRgErr=0.11,pH=7.0,fasta=list(fasta_P8D),ionic=0.15)
    proteins.loc['M9FP6Y'] = dict(eps_factor=0.2,temp=296,expRg=2.787,expRgErr=0.11,pH=7.0,fasta=list(fasta_M9FP6Y),ionic=0.15)
    proteins.loc['M9FP3Y'] = dict(eps_factor=0.2,temp=296,expRg=2.787,expRgErr=0.11,pH=7.0,fasta=list(fasta_M9FP3Y),ionic=0.15)
    proteins.loc['M12FP12Y'] = dict(eps_factor=0.2,temp=296,expRg=2.787,expRgErr=0.11,pH=7.0,fasta=list(fasta_M12FP12Y),ionic=0.15)
    proteins.loc['M8FP4Y'] = dict(eps_factor=0.2,temp=296,expRg=2.787,expRgErr=0.11,pH=7.0,fasta=list(fasta_M8FP4Y),ionic=0.15)
    proteins.loc['P7FM7Y'] = dict(eps_factor=0.2,temp=296,expRg=2.787,expRgErr=0.11,pH=7.0,fasta=list(fasta_P7FM7Y),ionic=0.15)
    proteins.loc['P7R'] = dict(eps_factor=0.2,temp=296,expRg=2.706,expRgErr=0.11,pH=7.0,fasta=list(fasta_P7R),ionic=0.15)
    proteins.loc['P2R'] = dict(eps_factor=0.2,temp=296,expRg=2.706,expRgErr=0.11,pH=7.0,fasta=list(fasta_P2R),ionic=0.15)
    proteins.loc['P7KP12D'] = dict(eps_factor=0.2,temp=296,expRg=2.706,expRgErr=0.11,pH=7.0,fasta=list(fasta_P7KP12D),ionic=0.15)
    proteins.loc['M3RP3K'] = dict(eps_factor=0.2,temp=296,expRg=2.706,expRgErr=0.11,pH=7.0,fasta=list(fasta_M3RP3K),ionic=0.15)
    proteins.loc['Ddx4WT'] = dict(eps_factor=0.2,temp=303,expRg=0,expRgErr=0,pH=6.5,fasta=list(fasta_Ddx4WT),ionic=0.13)
    proteins.loc['Ddx4CS'] = dict(eps_factor=0.2,temp=303,expRg=0,expRgErr=0,pH=6.5,fasta=list(fasta_Ddx4CS),ionic=0.13)
    proteins.loc['Ddx4RK'] = dict(eps_factor=0.2,temp=303,expRg=0,expRgErr=0,pH=6.5,fasta=list(fasta_Ddx4RK),ionic=0.13)
    proteins.loc['Ddx4FA'] = dict(eps_factor=0.2,temp=303,expRg=0,expRgErr=0,pH=6.5,fasta=list(fasta_Ddx4FA),ionic=0.13)
    proteins.loc['A2'] = dict(eps_factor=0.2,temp=298,expRg=0,expRgErr=0,pH=5.5,fasta=list(fasta_A2WT),ionic=0.010)
    proteins.loc['ht40'] = dict(eps_factor=0.2,temp=298,expRg=0,expRgErr=0,pH=7.4,fasta=list(fasta_ht40),ionic=0.005)
    proteins.loc['A2NS'] = dict(eps_factor=0.2,temp=298,expRg=0,expRgErr=0,pH=5.5,fasta=list(fasta_A2NS),ionic=0.010)
    return proteins


def genParamsLJ(df,name,prot):
    fasta = prot.fasta.copy()
    r = df.copy()
    r.loc['X'] = r.loc[fasta[0]]
    r.loc['Z'] = r.loc[fasta[-1]]
    r.loc['X','MW'] += 2
    r.loc['Z','MW'] += 16
    fasta[0] = 'X'
    fasta[-1] = 'Z'
    types = list(np.unique(fasta))
    MWs = [r.loc[a,'MW'] for a in types]
    sigmamap = pd.DataFrame((r.sigmas.values+r.sigmas.values.reshape(-1,1))/2,
                            index=r.sigmas.index,columns=r.sigmas.index)
    lambdamap = pd.DataFrame((r.lambdas.values+r.lambdas.values.reshape(-1,1))/2,
                             index=r.lambdas.index,columns=r.lambdas.index)
    lj_eps = prot.eps_factor*4.184
    # Generate pairs of amino acid types
    pairs = np.array(list(itertools.combinations_with_replacement(types,2)))
    return pairs, lj_eps, lambdamap, sigmamap, fasta, types, MWs

def genParamsDH(df,name,prot,temp):
    kT = 8.3145*temp*1e-3
    fasta = prot.fasta.copy()
    r = df.copy()
    # Set the charge on HIS based on the pH of the protein solution
    r.loc['H','q'] = 1. / ( 1 + 10**(prot.pH-6) )
    r.loc['X'] = r.loc[fasta[0]]
    r.loc['Z'] = r.loc[fasta[-1]]
    fasta[0] = 'X'
    fasta[-1] = 'Z'
    r.loc['X','q'] = r.loc[prot.fasta[0],'q'] + 1.
    r.loc['Z','q'] = r.loc[prot.fasta[-1],'q'] - 1.
    # Calculate the prefactor for the Yukawa potential
    qq = pd.DataFrame(r.q.values*r.q.values.reshape(-1,1),index=r.q.index,columns=r.q.index)
    fepsw = lambda T : 5321/T+233.76-0.9297*T+0.1417*1e-2*T*T-0.8292*1e-6*T**3
    epsw = fepsw(temp)
    lB = 1.6021766**2/(4*np.pi*8.854188*epsw)*6.022*1000/kT
    charges = [r.loc[a].q*np.sqrt(lB*kT) for a in fasta]
    yukawa_eps = qq*lB*kT
    # Calculate the inverse of the Debye length
    yukawa_kappa = np.sqrt(8*np.pi*lB*prot.ionic*6.022/10)
    return yukawa_eps, yukawa_kappa, charges

def genDCD(residues,name,prot,temp,n_chains):
    top = md.Topology()
    for _ in range(n_chains):
        chain = top.add_chain()
        for resname in prot.fasta:
            residue = top.add_residue(residues.loc[resname,'three'], chain)
            top.add_atom(residues.loc[resname,'three'], 
                         element=md.element.carbon, residue=residue)
        for i in range(chain.n_atoms-1):
            top.add_bond(chain.atom(i),chain.atom(i+1))

    t = md.load(name+'/{:d}'.format(temp)+'/{:s}.gsd'.format(name), top)
    lz = t.unitcell_lengths[0,2]
    edges = np.arange(-lz/2.,lz/2.,1)
    dz = (edges[1]-edges[0])/2.
    z = edges[:-1]+dz
    h = np.apply_along_axis(lambda a: np.histogram(a,bins=edges)[0], 1, t.xyz[:,:,2])
    #np.savetxt('{:s}_{:d}_na.dat'.format(name,temp),h)
    zmid = np.apply_along_axis(lambda a: z[a>np.quantile(a,.98)].mean(), 1, h)
    indices = np.argmin(np.abs(t.xyz[:,:,2]-zmid[:,np.newaxis]),axis=1)
    t[0].save_pdb(name+"/{:d}".format(temp)+'/top.pdb')
    t.save_dcd(name+"/{:d}".format(temp)+'/traj2.dcd')

    u = MDAnalysis.Universe(name+"/{:d}".format(temp)+'/top.pdb',name+'/{:d}'.format(temp)+'/traj2.dcd')
    ag = u.atoms
    with MDAnalysis.Writer(name+'/{:d}'.format(temp)+'/traj1.dcd', ag.n_atoms) as W:
        for ts,ndx in zip(u.trajectory,indices): 
            ts = transformations.unwrap(ag)(ts)
            ts = transformations.center_in_box(
                u.select_atoms('index {:d}'.format(ndx)), center='geometry')(ts)
            ts = transformations.wrap(ag)(ts)
            W.write(ag)
    
    t = md.load(name+'/{:d}'.format(temp)+'/traj1.dcd',top=name+'/{:d}'.format(temp)+'/top.pdb')
    lz = t.unitcell_lengths[0,2]
    edges = np.arange(0,lz,1)
    dz = (edges[1]-edges[0])/2.
    z = edges[:-1]+dz
    h = np.apply_along_axis(lambda a: np.histogram(a,bins=edges)[0], 1, t.xyz[:,:,2])
    h1 = np.mean(h[500:],axis=0)
    maxoverlap = np.apply_along_axis(lambda a: np.correlate(h1,np.histogram(a,
                bins=edges)[0], 'full').argmax()-h1.size+dz, 1, t.xyz[:,:,2])

    u = MDAnalysis.Universe(name+'/{:d}'.format(temp)+'/top.pdb',name+'/{:d}'.format(temp)+'/traj1.dcd')
    ag = u.atoms
    with MDAnalysis.Writer(name+'/{:d}'.format(temp)+'/traj.dcd', ag.n_atoms) as W:
        for ts,mo in zip(u.trajectory,maxoverlap):
            ts = transformations.unwrap(ag)(ts)
            ts = transformations.translate([0,0,mo*10])(ts)  
            ts = transformations.wrap(ag)(ts)
            W.write(ag)
   
    t = md.load(name+'/{:d}'.format(temp)+'/traj.dcd',top=name+'/{:d}'.format(temp)+'/top.pdb')
    h = np.apply_along_axis(lambda a: np.histogram(a,bins=edges)[0], 1, t.xyz[:,:,2])
    np.save('{:s}_{:d}.npy'.format(name,temp),h,allow_pickle=False)
    os.remove(name+'/{:d}'.format(temp)+'/traj1.dcd')
    os.remove(name+'/{:d}'.format(temp)+'/traj2.dcd')


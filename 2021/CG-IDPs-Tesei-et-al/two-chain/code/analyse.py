import pandas as pd
import numpy as np
import mdtraj as md
import itertools

def initProteinsDimers():
    proteins = pd.DataFrame(index=['FUS','FUS12E','A2','aSyn','ht40','p15PAF'], columns=['labels','eps_factor','wh','L','temp','obs','pH','ionic','expPREs','fasta','path'])
    fasta_FUS = """MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQ
SQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSS
SYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNS""".replace('\n', '')
    fasta_FUS12E = """GMASNDYEQQAEQSYGAYPEQPGQGYEQQSEQPYGQQSYSGYEQSTDTSGYGQSSYSSYGQ
EQNTGYGEQSTPQGYGSTGGYGSEQSEQSSYGQQSSYPGYGQQPAPSSTSGSYGSSEQSS
SYGQPQSGSYEQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNS""".replace('\n', '')
    fasta_A2 = """GHMGRGGNFGFGDSRGGGGNFGPGPGSNFRGGSDGYGSGRGFGDGYNGYGGGPGG
GNFGGSPGYGGGRGGYGGGGPGYGNQGGGYGGGYDNYGGGNYGSGNYNDFGNYNQQPSNYGPMKSGNFGGSRNMGG
PYGGGNYGPGGSGGSGGYGGRSRY""".replace('\n', '')
    fasta_aSyn = """MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTK
EQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA""".replace('\n', '') 
    fasta_ht40 = """MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEP
GSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKK
AKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREP
KKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGS
VQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTD
HGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL""".replace('\n', '') 
    fasta_p15PAF = """MVRTKADSVPGTYRKVVAARAPRKVLGSSTSATNSTSVSSRKAENKYAGGNPVCVRPTPK
WQKGIGEFFRLSPKDSEKENQIPEEAGSSGLGKAKRKACPLQPDHTNDEKE""".replace('\n', '') 
    proteins.loc['FUS'] = dict(labels=[16, 86, 142],eps_factor=0.2,L=40.5,wh=850,temp=298,obs='rate',pH=5.5,fasta=list(fasta_FUS),ionic=0.15,path='ff') 
    proteins.loc['FUS12E'] = dict(labels=[16, 86, 142],eps_factor=0.2,L=40.5,wh=850,temp=298,obs='rate',pH=5.5,fasta=list(fasta_FUS12E),ionic=0.15,path='ff')
    proteins.loc['A2'] = dict(labels=[99, 143],eps_factor=0.2,L=48,wh=850,temp=298,obs='rate',pH=5.5,fasta=list(fasta_A2),ionic=0.005,path='ff')
    proteins.loc['aSyn'] = dict(eps_factor=0.2,temp=283,pH=7.4,ionic=0.125,L=25.5,fasta=list(fasta_aSyn),path='ff')
    proteins.loc['ht40'] = dict(eps_factor=0.2,temp=278,pH=6.8,ionic=0.100,L=48.0,fasta=list(fasta_ht40),path='ff')
    proteins.loc['p15PAF'] = dict(eps_factor=0.2,temp=298,pH=7.0,ionic=0.150,L=34.0,fasta=list(fasta_p15PAF),path='ff')
    return proteins

def genParamsLJ(df,name,prot):
    fasta = prot.fasta.copy()
    r = df.copy()
    r.loc['X'] = r.loc[fasta[0]]
    r.loc['Z'] = r.loc[fasta[-1]]
    fasta[0] = 'X'
    fasta[-1] = 'Z'
    types = list(np.unique(fasta))
    sigmamap = pd.DataFrame((r.sigmas.values+r.sigmas.values.reshape(-1,1))/2,
                            index=r.sigmas.index,columns=r.sigmas.index)
    lambdamap = pd.DataFrame((r.lambdas.values+r.lambdas.values.reshape(-1,1))/2,
                             index=r.lambdas.index,columns=r.lambdas.index)
    lj_eps = prot.eps_factor*4.184
    # Generate pairs of amino acid types
    pairs = np.array(list(itertools.combinations_with_replacement(types,2)))
    return pairs, lj_eps, lambdamap, sigmamap, fasta, types

def genParamsDH(df,name,prot):
    kT = 8.3145*prot.temp*1e-3
    r = df.copy()
    # Set the charge on HIS based on the pH of the protein solution
    r.loc['H','q'] = 1. / ( 1 + 10**(prot.pH-6) )
    r.loc['X','q'] = r.loc[prot.fasta[0],'q'] + 1.
    r.loc['Z','q'] = r.loc[prot.fasta[-1],'q'] - 1.
    # Calculate the prefactor for the Yukawa potential
    qq = pd.DataFrame(r.q.values*r.q.values.reshape(-1,1),index=r.q.index,columns=r.q.index)
    fepsw = lambda T : 5321/T+233.76-0.9297*T+0.1417*1e-2*T*T-0.8292*1e-6*T**3
    epsw = fepsw(prot.temp)
    lB = 1.6021766**2/(4*np.pi*8.854188*epsw)*6.022*1000/kT
    yukawa_eps = qq*lB*kT
    # Calculate the inverse of the Debye length
    yukawa_kappa = np.sqrt(8*np.pi*lB*prot.ionic*6.022/10)
    return yukawa_eps, yukawa_kappa

def genDCD(residues,name,prot,path,run_type,n_chains):
    """ 
    Generates coordinate and trajectory 
    in convenient formats for multiple chains
    """
    top = md.Topology()
    for _ in range(n_chains):
        chain = top.add_chain()
        for resname in prot.fasta:
            residue = top.add_residue(residues.loc[resname,'three'], chain)
            top.add_atom(residues.loc[resname,'three'], element=md.element.carbon, residue=residue)
    traj = md.load_dcd(path+"/{:s}.dcd".format(run_type), top)
    traj.center_coordinates()
    traj.xyz *= 10
    traj.unitcell_lengths *= 10
    traj.xyz += traj.unitcell_lengths[0,0]/2
    traj[:].save_dcd(path+"/{:s}.dcd".format(name))
    traj[0].save_pdb(path+"/{:s}.pdb".format(name))

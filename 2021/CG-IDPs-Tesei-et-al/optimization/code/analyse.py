import pandas as pd
import numpy as np
import mdtraj as md
import itertools
from DEERPREdict.utils import Operations

def loadExpPREs(name,prot):
    value = {}
    error = {}
    resnums = np.arange(1,len(prot.fasta)+1)
    for label in prot.labels:
        value[label], error[label] = np.loadtxt(name+'/expPREs/exp-{:d}.dat'.format(label),unpack=True)
    v = pd.DataFrame(value,index=resnums)
    v.rename_axis('residue', inplace=True)
    v.rename_axis('label', axis='columns',inplace=True)
    e = pd.DataFrame(error,index=resnums)
    e.rename_axis('residue', inplace=True)
    e.rename_axis('label', axis='columns',inplace=True)
    return pd.concat(dict(value=v,error=e),axis=1)

def loadInitPREs(name,prot):
    obs = 1 if prot.obs=='ratio' else 2
    value = {}
    resnums = np.arange(1,len(prot.fasta)+1)
    for label in prot.labels:
        value[label] = np.loadtxt(prot.path+'/calcPREs/res-{:d}.dat'.format(label))[:,obs]
    v = pd.DataFrame(value,index=resnums)
    v.rename_axis('residue', inplace=True)
    v.rename_axis('label', axis='columns',inplace=True)
    return v

def calcChi2(prot):
    obs = 1 if prot.obs=='ratio' else 2
    chi2 = 0
    for label in prot.labels:
        y = np.loadtxt(prot.path+'/calcPREs/res-{:d}.dat'.format(label))[:,obs]
        chi = (prot.expPREs.value[label].values - y) / prot.expPREs.error[label].values
        chi = chi[~np.isnan(chi)]
        chi2 += np.nansum( np.power( chi, 2) ) / chi.size
    return 0.1*chi2 / len(prot.labels)

def optTauC(prot):
    obs = 1 if prot.obs=='ratio' else 2
    chi2list = []
    tau_c = np.arange(1,15.05,1)
    for tc in tau_c:
        chi2 = 0
        for label in prot.labels:
            x,y = np.loadtxt(prot.path+'/calcPREs/res-{:d}.dat'.format(label),usecols=(0,1),unpack=True)
            measured_resnums = np.where(~np.isnan(y))[0] 
            data = pd.read_pickle(prot.path+'/calcPREs/res-{:d}.pkl'.format(label), compression='gzip')
            gamma_2_av = np.full(y.size, fill_value=np.NaN)
            s_pre = np.power(data['r3'], 2)/data['r6']*data['angular']
            gamma_2 = Operations.calc_gamma_2(data['r6'], s_pre, tau_c = tc * 1e-9, tau_t = 1e-10, wh = prot.wh, k = 1.23e16)
            gamma_2 = np.ma.MaskedArray(gamma_2, mask = np.isnan(gamma_2))
            gamma_2_av[measured_resnums] = np.ma.average(gamma_2, axis=0).data
            if prot.obs == 'ratio':
                y = 10 * np.exp(-gamma_2_av * 0.01) / ( 10 + gamma_2_av )
            else:
                y = gamma_2_av
            chi = (prot.expPREs.value[label].values - y) / prot.expPREs.error[label].values
            chi = chi[~np.isnan(chi)]
            chi2 += np.nansum( np.power( chi, 2) ) / chi.size
        chi2list.append(0.1 * chi2 / len(prot.labels))

    tc_min = tau_c[np.argmin(chi2list)]

    for label in prot.labels:
        x,y = np.loadtxt(prot.path+'/calcPREs/res-{:d}.dat'.format(label),usecols=(0,1),unpack=True)
        measured_resnums = np.where(~np.isnan(y))[0] 
        data = pd.read_pickle(prot.path+'/calcPREs/res-{:d}.pkl'.format(label), compression='gzip')
        gamma_2_av = np.full(y.size, fill_value=np.NaN)
        s_pre = np.power(data['r3'], 2)/data['r6']*data['angular']
        gamma_2 = Operations.calc_gamma_2(data['r6'], s_pre, tau_c = tc_min * 1e-9, tau_t = 1e-10, wh = prot.wh, k = 1.23e16)
        gamma_2 = np.ma.MaskedArray(gamma_2, mask = np.isnan(gamma_2))
        gamma_2_av[measured_resnums] = np.ma.average(gamma_2, axis=0).data
        i_ratio = 10 * np.exp(-gamma_2_av * 0.01) / ( 10 + gamma_2_av )
        np.savetxt(prot.path+'/calcPREs/res-{:d}.dat'.format(label),np.c_[x,i_ratio,gamma_2_av])
 
    return tc_min, calcChi2(prot)

def Rg2invRh(x,N,a1=0.216,a2=4.06,a3=0.821):
    # convert Rg to Rh (DOI: 10.1016/j.bpj.2017.06.042)
    return (a1*(x*10-a2*N**0.33)/(N**0.6-N**0.33)+a3)/x

def calcRh(df,name,prot):
    if isinstance(prot.rgarray, np.ndarray):
        rg = np.sqrt( np.dot(np.power(prot.rgarray,2),prot.weights) )
        rh = 1/np.dot(Rg2invRh(prot.rgarray,len(prot.fasta)),prot.weights)
        rhkr = 1/np.dot(prot.invrij,prot.weights)
        chi2_rh = np.power((prot.expRh-rhkr)/(rh-rhkr)*2,2)
        return rg, rh, rhkr, chi2_rh
    else:
        t = md.load_dcd(prot.path+"/{:s}.dcd".format(name),prot.path+"/{:s}.pdb".format(name))
        residues = [res.name for res in t.top.atoms]
        masses = df.loc[residues,'MW'].values
        masses[0] += 2
        masses[-1] += 16
        rgarray = md.compute_rg(t,masses=masses)
        invrij = (1-1/len(residues))*(1/md.compute_distances(t,t.top.select_pairs('all','all'))).mean(axis=1)
        rg = np.sqrt( np.power(rgarray,2).mean() )
        rh = 1/Rg2invRh(rgarray,len(residues)).mean()
        rhkr = 1/invrij.mean()
        chi2_rh = np.power((prot.expRh-rhkr)/(rh-rhkr)*2,2)
        return rgarray, invrij, rg, rh, rhkr, chi2_rh

def calcRg(df,name,prot):
    if isinstance(prot.rgarray, np.ndarray):
        rg = np.sqrt( np.dot(np.power(prot.rgarray,2), prot.weights) )
        chi2_rg = np.power((prot.expRg-rg)/prot.expRgErr,2)
        return rg, chi2_rg
    else:
        t = md.load_dcd(prot.path+"/{:s}.dcd".format(name),prot.path+"/{:s}.pdb".format(name))
        residues = [res.name for res in t.top.atoms]
        masses = df.loc[residues,'MW'].values
        rgarray = md.compute_rg(t,masses=masses)
        rg = np.sqrt( np.power(rgarray, 2).mean() )
        chi2_rg = np.power((prot.expRg-rg)/prot.expRgErr,2)
        return rgarray, rg, chi2_rg

def initProteins(outdir):
    outdir = '/'+outdir if len(outdir)>0 else outdir
    proteins = pd.DataFrame(index=['OPN','FUS','FUS12E','Sic1','aSyn','A2'],
            columns=['labels','eps_factor','wh','tau_c','temp','expRh','Rh','RhKR','Rg','obs','pH','rgarray','ionic','invrij','expPREs','initPREs','eff','chi2_pre','chi2_rh','fasta','weights','path'])
    fasta_OPN = """MHQDHVDSQSQEHLQQTQNDLASLQQTHYSSEENADVPEQPDFPDV
PSKSQETVDDDDDDDNDSNDTDESDEVFTDFPTEAPVAPFNRGDNAGRGDSVAYGFRAKA
HVVKASKIRKAARKLIEDDATTEDGDSQPAGLWWPKESREQNSRELPQHQSVENDSRPKF
DSREVDGGDSKASAGVDSRESQGSVPAVDASNQTLESAEDAEDRHSIENNEVTR""".replace('\n', '')
    fasta_FUS = """MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQ
SQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSS
SYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNS""".replace('\n', '')
    fasta_FUS12E = """GMASNDYEQQAEQSYGAYPEQPGQGYEQQSEQPYGQQSYSGYEQSTDTSGYGQSSYSSYGQ
EQNTGYGEQSTPQGYGSTGGYGSEQSEQSSYGQQSSYPGYGQQPAPSSTSGSYGSSEQSS
SYGQPQSGSYEQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNS""".replace('\n', '')
    fasta_Sic1 = """GSMTPSTPPRSRGTRYLAQPSGNTSSSALMQGQKTPQKPS
QNLVPVTPSTTKSFKNAPLLAPPNSNMGMTSPFNGLTSPQRSPFPKSSVKRT""".replace('\n', '')
    fasta_aSyn = """MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTK
EQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDP
DNEAYEMPSEEGYQDYEPEA""".replace('\n', '')
    fasta_A2 = """GHMGRGGNFGFGDSRGGGGNFGPGPGSNFRGGSDGYGSGR
GFGDGYNGYGGGPGGGNFGGSPGYGGGRGGYGGGGPGYGNQGGGYGGGYDNYGGGNYGSG
NYNDFGNYNQQPSNYGPMKSGNFGGSRNMGGPYGGGNYGPGGSGGSGGYGGRSRY""".replace('\n', '')

    proteins.loc['A2'] = dict(labels=[99, 143], eps_factor=0.2,tau_c=10.0,
                               wh=850,temp=298,expRh=2.89,obs='ratio',pH=5.5,fasta=list(fasta_A2),ionic=0.005,weights=False,path='A2'+outdir)
    proteins.loc['aSyn'] = dict(labels=[24, 42, 62, 87, 103], eps_factor=0.2,tau_c=1.0,
                               wh=700,temp=283,expRh=2.82,obs='ratio',pH=7.4,fasta=list(fasta_aSyn),ionic=0.15,weights=False,path='aSyn'+outdir)
    proteins.loc['OPN'] = dict(labels=[10, 33, 64, 88, 117, 130, 144, 162, 184, 203], eps_factor=0.2,tau_c=3.0,
                               wh=800,temp=298,expRh=4.107,obs='rate',pH=6.5,fasta=list(fasta_OPN),ionic=0.15,weights=False,path='OPN'+outdir)
    proteins.loc['FUS'] = dict(labels=[16, 86, 142], eps_factor=0.2,tau_c=10.0,
                               wh=850,temp=298,expRh=3.32,obs='rate',pH=5.5,fasta=list(fasta_FUS),ionic=0.15,weights=False,path='FUS'+outdir)
    proteins.loc['FUS12E'] = dict(labels=[16, 86, 142], eps_factor=0.2,tau_c=10.0,
                               wh=850,temp=298,expRh=3.97,obs='rate',pH=5.5,fasta=list(fasta_FUS12E),ionic=0.15,weights=False,path='FUS12E'+outdir)
    proteins.loc['Sic1'] = dict(labels=[1, 23, 40, 66, 85, 92], eps_factor=0.2,tau_c=8.0,
                               wh=500,temp=278,expRh=2.15,obs='rate',pH=7.0,fasta=list(fasta_Sic1),ionic=0.15,weights=False,path='Sic1'+outdir)
    return proteins

def initProteinsRgs(outdir):
    outdir = '/'+outdir if len(outdir)>0 else outdir
    proteins = pd.DataFrame(index=['A1','SH4UD','Hst5','aSyn140','PNt','Hst52','ACTR','RNaseA','p15PAF','OPN220','Sic92','FhuA','CoRNID','ColNT','hNL3cyt','K10','K27','K25','K32','K23','K44','M12FP12Y','P7FM7Y','M9FP6Y','M8FP4Y','M9FP3Y','M10R','M6R','P2R','P7R','M3RP3K','M6RP6K','M10RP10K','M4D','P4D','P8D','P12D','P12E','P7KP12D','P7KP12Db','M12FP12YM10R','M10FP7RP12D'],
            columns=['eps_factor','temp','expRg','expRgErr','Rg','rgarray','eff','chi2_rg','weights','pH','ionic','fasta','path'])
    fasta_Hst52 = """DSHAKRHHGYKRKFHEKHHSHRGYDSHAKRHHGYKRKFHEKHHSHRGY""".replace('\n', '') # DOI: 10.1021/acs.jpcb.0c09635
    fasta_Hst5 = """DSHAKRHHGYKRKFHEKHHSHRGY""".replace('\n', '') 
    fasta_aSyn140 = """MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTK
EQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA""".replace('\n', '')
    fasta_PNt = """DWNNQSIVKTGERQHGIHIQGSDPGGVRTASGTTIKVSGRQAQGILLENPAAELQFRNGSVTSSGQLSDDGIRRFLGTVTVKAGKLVADHATLANVGDTWDDDGIALYVAGEQAQASIAD
STLQGAGGVQIERGANVTVQRSAIVDGGLHIGALQSLQPEDLPPSRVVLRDTNVTAVPASGAPAAVSVLGASELTLDGGHITGGRAAGVAAMQGAVVHLQRATIRRGDALAGGAVPGGAVPGGAVPGGFGPG
GFGPVLDGWYGVDVSGSSVELAQSIVEAPELGAAIRVGRGARVTVPGGSLSAPHGNVIETGGARRFAPQAAPLSITLQAGAH""".replace('\n', '')
    fasta_ACTR = """GTQNRPLLRNSLDDLVGPPSNLEGQSDERALLDQLHTLLSNTDATGLEEIDRALGIPELVNQGQALEPKQD""".replace('\n', '') # DOI: 10.1073/pnas.1322611111
    fasta_RNaseA = """KETAAAKFERQHMDSSTSAASSSNYCNQMMKSRNLTKDRCKPVNTFVHESLADVQAVCSQKNVACKNGQTNCYQSYSTMSI
TDCRETGSSKYPNCAYKTTQANKHIIVACEGNPYVPVHFDASV""".replace('\n', '')
    fasta_p15PAF = """MVRTKADSVPGTYRKVVAARAPRKVLGSSTSATNSTSVSSRKAENKYAGGNPVCVRPTPKWQKGIGEFFR
LSPKDSEKENQIPEEAGSSGLGKAKRKACPLQPDHTNDEKE""".replace('\n', '') # DOI: 10.1016/j.bpj.2013.12.046
    fasta_CoRNID = """GPHMQVPRTHRLITLADHICQIITQDFARNQVPSQASTSTFQTSPSALSSTPVRTKTSSRYS
PESQSQTVLHPRPGPRVSPENLVDKSRGSRPGKSPERSHIPSEPYEPISPPQGPAVHEKQDSMLLLSQRGVDPAEQRSDSRSP
GSISYLPSFFTKLESTSPMVKSKKQEIFRKLNSSGGGDSDMAAAQPGTEIFNLPAVTTSGAVSSRSHSFADPASNLGLEDIIR
KALMGSFDDKVEDHGVVMSHPVGIMPGSASTSVVTSSEARRDE""".replace('\n', '') # SASDF34
    fasta_OPN220 = """MHQDHVDSQSQEHLQQTQNDLASLQQTHYSSEENADVPEQPDFPDV
PSKSQETVDDDDDDDNDSNDTDESDEVFTDFPTEAPVAPFNRGDNAGRGDSVAYGFRAKA
HVVKASKIRKAARKLIEDDATTEDGDSQPAGLWWPKESREQNSRELPQHQSVENDSRPKF
DSREVDGGDSKASAGVDSRESQGSVPAVDASNQTLESAEDAEDRHSIENNEVTR""".replace('\n', '')
    fasta_Sic92 = """GSMTPSTPPRSRGTRYLAQPSGNTSSSALMQGQKTPQKPS
QNLVPVTPSTTKSFKNAPLLAPPNSNMGMTSPFNGLTSPQRSPFPKSSVKRT""".replace('\n', '')
    fasta_FhuA = """SESAWGPAATIAARQSATGTKTDTPIQKVPQSISVVTAEEMALHQPKSVKEALSYTPGVSVGTRGASNTYDHLIIRGFAAEGQS
QNNYLNGLKLQGNFYNDAVIDPYMLERAEIMRGPVSVLYGKSSPGGLLNMVSKRPTTEPL""".replace('\n', '')
    fasta_K44 =  """MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKA
KGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRL
QTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIE""".replace('\n', '')
    fasta_K10 = """MQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVK
SEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL""".replace('\n', '')
    fasta_K27 = """MSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIV
YKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVY""".replace('\n', '')
    fasta_K25 = """MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDK
KAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRL""".replace('\n', '')
    fasta_K32 = """MSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQII
NKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVY""".replace('\n', '')
    fasta_K23 = """MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDK
KAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRL
THKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL""".replace('\n', '')
    fasta_A1 = """GSMASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFGGRSSGGSGGGGQYFAKPRNQGGYGGSSSSSSYGSGRRF""".replace('\n', '')
    fasta_M12FP12Y = """GSMASASSSQRGRSGSGNYGGGRGGGYGGNDNYGRGGNYSGRGGYGGSRGGGGYGGSGDGYNGYGNDGSNYGGGGSYNDYGNYNNQ
SSNYGPMKGGNYGGRSSGGSGGGGQYYAKPRNQGGYGGSSSSSSYGSGRRY""".replace('\n', '')
    fasta_P7FM7Y = """GSMASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGFGGSGDGFNGFGNDGSNFGGGGSFNDFGNFNNQ
SSNFGPMKGGNFGGRSSGGSGGGGQFFAKPRNQGGFGGSSSSSSFGSGRRF""".replace('\n', '')
    fasta_M9FP6Y = """GSMASASSSQRGRSGSGNFGGGRGGGYGGNDNYGRGGNYSGRGGFGGSRGGGGYGGSGDGYNGGGNDGSNYGGGGSYNDSGNYNNQ
SSNFGPMKGGNYGGRSSGGSGGGGQYGAKPRNQGGYGGSSSSSSYGSGRRY""".replace('\n', '')
    fasta_M8FP4Y = """GSMASASSSQRGRSGSGNFGGGRGGGYGGNDNGGRGGNYSGRGGFGGSRGGGGYGGSGDGYNGGGNDGSNYGGGGSYNDSGNYNNQ
SSNFGPMKGGNYGGRSSGGSGGGGQYGAKPRNQGGYGGSSSSSSYGSGRRF""".replace('\n', '')
    fasta_M9FP3Y = """GSMASASSSQRGRSGSGNFGGGRGGGYGGNDNGGRGGNYSGRGGFGGSRGGGGYGGSGDGYNGGGNDGSNYGGGGSYNDSGNGNNQ
SSNFGPMKGGNYGGRSSGGSGGGGQYGAKPRNQGGYGGSSSSSSYGSGRRS""".replace('\n', '')
    fasta_M10R = """GSMASASSSQGGSSGSGNFGGGGGGGFGGNDNFGGGGNFSGSGGFGGSGGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFGGSSSGPYGGGGQYFAKPGNQGGYGGSSSSSSYGSGGGF""".replace('\n', '')
    fasta_M6R = """GSMASASSSQGGRSGSGNFGGGRGGGFGGNDNFGGGGNFSGSGGFGGSRGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFGGSSSGPYGGGGQYFAKPGNQGGYGGSSSSSSYGSGGRF""".replace('\n', '')
    fasta_P2R = """GSMASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGYGGSGDGYNGFRNDGSNFGGGGRYNDFGNYNNQ
SSNFGPMKGGNFGGRSSGPYGGGGQYFAKPRNQGGYGGSSSSSSYGSGRRF""".replace('\n', '')
    fasta_P7R = """GSMASASSSQRGRSGRGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGRYGGSGDRYNGFGNDGRNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFRGRSSGPYGRGGQYFAKPRNQGGYGGSSSSRSYGSGRRF""".replace('\n', '')
    fasta_M3RP3K = """GSMASASSSQRGKSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSKGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFGGRSSGGSGGGGQYFAKPRNQGGYGGSSSSSSYGSGRKF""".replace('\n', '')
    fasta_M6RP6K = """GSMASASSSQKGKSGSGNFGGGRGGGFGGNDNFGKGGNFSGRGGFGGSKGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFGGKSSGGSGGGGQYFAKPRNQGGYGGSSSSSSYGSGRKF""".replace('\n', '')
    fasta_M10RP10K = """GSMASASSSQKGKSGSGNFGGGKGGGFGGNDNFGKGGNFSGKGGFGGSKGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFGGKSSGGSGGGGQYFAKPKNQGGYGGSSSSSSYGSGKKF""".replace('\n', '')
    fasta_M4D = """GSMASASSSQRGRSGSGNFGGGRGGGFGGNGNFGRGGNFSGRGGFGGSRGGGGYGGSGGGYNGFGNSGSNFGGGGSYNGFGNYNNQ
SSNFGPMKGGNFGGRSSGPYGGGGQYFAKPRNQGGYGGSSSSSSYGSGRRF""".replace('\n', '')
    fasta_P4D = """GSMASASSSQRDRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGDFGGSRGGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFGGRSSDPYGGGGQYFAKPRNQGGYGGSSSSSSYDSGRRF""".replace('\n', '')
    fasta_P8D = """GSMASASSSQRDRSGSGNFGGGRDGGFGGNDNFGRGDNFSGRGDFGGSRDGGGYGGSGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFGPMKGGNFGGRSSDPYGGGGQYFAKPRNQDGYGGSSSSSSYDSGRRF""".replace('\n', '')
    fasta_P12D = """GSMASADSSQRDRDDSGNFGDGRGGGFGGNDNFGRGGNFSDRGGFGGSRGDGGYGGDGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFDPMKGGNFGDRSSGPYDGGGQYFAKPRNQGGYGGSSSSSSYGSDRRF""".replace('\n', '')
    fasta_P12E = """GSMASAESSQREREESGNFGEGRGGGFGGNDNFGRGGNFSERGGFGGSRGEGGYGGEGDGYNGFGNDGSNFGGGGSYNDFGNYNNQ
SSNFEPMKGGNFGERSSGPYEGGGQYFAKPRNQGGYGGSSSSSSYGSERRF""".replace('\n', '')
    fasta_P7KP12D = """GSMASADSSQRDRDDKGNFGDGRGGGFGGNDNFGRGGNFSDRGGFGGSRGDGKYGGDGDKYNGFGNDGKNFGGGGSYNDFGNYNNQ
SSNFDPMKGGNFKDRSSGPYDKGGQYFAKPRNQGGYGGSSSSKSYGSDRRF""".replace('\n', '')
    fasta_P7KP12Db = """GSMASAKSSQRDRDDDGNFGKGRGGGFGGNKNFGRGGNFSKRGGFGGSRGKGKYGGKGDDYNGFGNDGDNFGGGGSYNDFGNYNNQ
SSNFDPMDGGNFDDRSSGPYDDGGQYFADPRNQGGYGGSSSSKSYGSKRRF""".replace('\n', '')
    fasta_M12FP12YM10R = """GSMASASSSQGGSSGSGNYGGGGGGGYGGNDNYGGGGNYSGSGGYGGSGGGGGYGGSGDGYNGYGNDGSNYGGGGSYNDYGNYNNQ
SSNYGPMKGGNYGGSSSGPYGGGGQYYAKPGNQGGYGGSSSSSSYGSGGGY""".replace('\n', '')
    fasta_M10FP7RP12D = """GSMASADSSQRDRDDRGNFGDGRGGGGGGNDNFGRGGNGSDRGGGGGSRGDGRYGGDGDRYNGGGNDGRNGGGGGSYNDGGNYNNQ
SSNGDPMKGGNGRDRSSGPYDRGGQYGAKPRNQGGYGGSSSSRSYGSDRRG""".replace('\n', '')
    fasta_SH4UD = """MGSNKSKPKDASQRRRSLEPAENVHGAGGGAFPASQTPSKPASADGHRGPSAAFAPAAA
EPKLFGGFNSSDTVTSPQRAGPLAGGSAWSHPQFEK""".replace('\n', '') # DOI: 
    fasta_hNL3cyt = """MYRKDKRRQEPLRQPSPQRGAWAPELGAAPEEELAALQLGPTHHECEAG
PPHDTLRLTALPDYTLTLRRSPDDIPLMTPNTITMIPNSLVG
LQTLHPYNTFAAGFNSTGLPHSHSTTRV""".replace('\n', '') # DOI: 10.1529/biophysj.107.126995
    fasta_ColNT = """MGSNGADNAHNNAFGGGKNPGIGNTSGAGSNGSASSNRGNSNGWSWSNKPHKNDGFHSDGSYHITFHGDNNSKPKPGGNSGNRGNNGDGASSHHHHHH""".replace('\n', '') # SASDC53

    proteins.loc['ColNT'] = dict(eps_factor=0.2,temp=277,expRg=2.83,expRgErr=0.1,pH=7.6,fasta=list(fasta_ColNT),ionic=0.4,path='ColNT'+outdir)
    proteins.loc['hNL3cyt'] = dict(eps_factor=0.2,temp=293,expRg=3.15,expRgErr=0.2,pH=8.5,fasta=list(fasta_hNL3cyt),ionic=0.3,path='hNL3cyt'+outdir)
    proteins.loc['SH4UD'] = dict(eps_factor=0.2,temp=293,expRg=2.71,expRgErr=0.1,pH=8.0,fasta=list(fasta_SH4UD),ionic=0.2,path='SH4UD'+outdir)
    proteins.loc['OPN220'] = dict(eps_factor=0.2,temp=298,expRg=5.13,expRgErr=0.2,pH=6.5,fasta=list(fasta_OPN220),ionic=0.15,path='OPN220'+outdir)
    proteins.loc['Sic92'] = dict(eps_factor=0.2,temp=293,expRg=3.0,expRgErr=0.4,pH=7.5,fasta=list(fasta_Sic92),ionic=0.2,path='Sic92'+outdir)
    proteins.loc['FhuA'] = dict(eps_factor=0.2,temp=298,expRg=3.34,expRgErr=0.1,pH=7.5,fasta=list(fasta_FhuA),ionic=0.15,path='FhuA'+outdir)
    proteins.loc['K10'] = dict(eps_factor=0.2,temp=288,expRg=4.0,expRgErr=0.1,pH=7.4,fasta=list(fasta_K10),ionic=0.15,path='K10'+outdir)
    proteins.loc['K27'] = dict(eps_factor=0.2,temp=288,expRg=3.7,expRgErr=0.2,pH=7.4,fasta=list(fasta_K27),ionic=0.15,path='K27'+outdir)
    proteins.loc['K25'] = dict(eps_factor=0.2,temp=288,expRg=4.1,expRgErr=0.2,pH=7.4,fasta=list(fasta_K25),ionic=0.15,path='K25'+outdir)
    proteins.loc['K32'] = dict(eps_factor=0.2,temp=288,expRg=4.2,expRgErr=0.3,pH=7.4,fasta=list(fasta_K32),ionic=0.15,path='K32'+outdir)
    proteins.loc['K23'] = dict(eps_factor=0.2,temp=288,expRg=4.9,expRgErr=0.2,pH=7.4,fasta=list(fasta_K23),ionic=0.15,path='K23'+outdir)
    proteins.loc['K44'] = dict(eps_factor=0.2,temp=288,expRg=5.2,expRgErr=0.2,pH=7.4,fasta=list(fasta_K44),ionic=0.15,path='K44'+outdir)
    proteins.loc['A1'] = dict(eps_factor=0.2,temp=298,expRg=2.760,expRgErr=0.05,pH=7.0,fasta=list(fasta_A1),ionic=0.15,path='A1'+outdir)
    proteins.loc['M12FP12Y'] = dict(eps_factor=0.2,temp=298,expRg=2.604,expRgErr=0.05,pH=7.0,fasta=list(fasta_M12FP12Y),ionic=0.15,path='M12FP12Y'+outdir)
    proteins.loc['P7FM7Y'] = dict(eps_factor=0.2,temp=298,expRg=2.718,expRgErr=0.05,pH=7.0,fasta=list(fasta_P7FM7Y),ionic=0.15,path='P7FM7Y'+outdir) 
    proteins.loc['M9FP6Y'] = dict(eps_factor=0.2,temp=298,expRg=2.655,expRgErr=0.05,pH=7.0,fasta=list(fasta_M9FP6Y),ionic=0.15,path='M9FP6Y'+outdir) 
    proteins.loc['M8FP4Y'] = dict(eps_factor=0.2,temp=298,expRg=2.707,expRgErr=0.05,pH=7.0,fasta=list(fasta_M8FP4Y),ionic=0.15,path='M8FP4Y'+outdir) 
    proteins.loc['M9FP3Y'] = dict(eps_factor=0.2,temp=298,expRg=2.683,expRgErr=0.05,pH=7.0,fasta=list(fasta_M9FP3Y),ionic=0.15,path='M9FP3Y'+outdir) 
    proteins.loc['M10R'] = dict(eps_factor=0.2,temp=298,expRg=2.671,expRgErr=0.05,pH=7.0,fasta=list(fasta_M10R),ionic=0.15,path='M10R'+outdir)
    proteins.loc['M6R'] = dict(eps_factor=0.2,temp=298,expRg=2.573,expRgErr=0.05,pH=7.0,fasta=list(fasta_M6R),ionic=0.15,path='M6R'+outdir) 
    proteins.loc['P2R'] = dict(eps_factor=0.2,temp=298,expRg=2.623,expRgErr=0.05,pH=7.0,fasta=list(fasta_P2R),ionic=0.15,path='P2R'+outdir) 
    proteins.loc['P7R'] = dict(eps_factor=0.2,temp=298,expRg=2.709,expRgErr=0.05,pH=7.0,fasta=list(fasta_P7R),ionic=0.15,path='P7R'+outdir) 
    proteins.loc['M3RP3K'] = dict(eps_factor=0.2,temp=298,expRg=2.634,expRgErr=0.05,pH=7.0,fasta=list(fasta_M3RP3K),ionic=0.15,path='M3RP3K'+outdir) 
    proteins.loc['M6RP6K'] = dict(eps_factor=0.2,temp=298,expRg=2.787,expRgErr=0.05,pH=7.0,fasta=list(fasta_M6RP6K),ionic=0.15,path='M6RP6K'+outdir) 
    proteins.loc['M10RP10K'] = dict(eps_factor=0.2,temp=298,expRg=2.849,expRgErr=0.05,pH=7.0,fasta=list(fasta_M10RP10K),ionic=0.15,path='M10RP10K'+outdir)
    proteins.loc['M4D'] = dict(eps_factor=0.2,temp=298,expRg=2.642,expRgErr=0.05,pH=7.0,fasta=list(fasta_M4D),ionic=0.15,path='M4D'+outdir) 
    proteins.loc['P4D'] = dict(eps_factor=0.2,temp=298,expRg=2.718,expRgErr=0.05,pH=7.0,fasta=list(fasta_P4D),ionic=0.15,path='P4D'+outdir) 
    proteins.loc['P8D'] = dict(eps_factor=0.2,temp=298,expRg=2.685,expRgErr=0.05,pH=7.0,fasta=list(fasta_P8D),ionic=0.15,path='P8D'+outdir) 
    proteins.loc['P12D'] = dict(eps_factor=0.2,temp=298,expRg=2.801,expRgErr=0.05,pH=7.0,fasta=list(fasta_P12D),ionic=0.15,path='P12D'+outdir)
    proteins.loc['P12E'] = dict(eps_factor=0.2,temp=298,expRg=2.852,expRgErr=0.05,pH=7.0,fasta=list(fasta_P12E),ionic=0.15,path='P12E'+outdir)
    proteins.loc['P7KP12D'] = dict(eps_factor=0.2,temp=298,expRg=2.921,expRgErr=0.05,pH=7.0,fasta=list(fasta_P7KP12D),ionic=0.15,path='P7KP12D'+outdir) 
    proteins.loc['P7KP12Db'] = dict(eps_factor=0.2,temp=298,expRg=2.562,expRgErr=0.05,pH=7.0,fasta=list(fasta_P7KP12Db),ionic=0.15,path='P7KP12Db'+outdir)
    proteins.loc['M12FP12YM10R'] = dict(eps_factor=0.2,temp=298,expRg=2.607,expRgErr=0.05,pH=7.0,fasta=list(fasta_M12FP12YM10R),ionic=0.15,path='M12FP12YM10R'+outdir)
    proteins.loc['M10FP7RP12D'] = dict(eps_factor=0.2,temp=298,expRg=2.860,expRgErr=0.05,pH=7.0,fasta=list(fasta_M10FP7RP12D),ionic=0.15,path='M10FP7RP12D'+outdir) 
    proteins.loc['Hst5'] = dict(eps_factor=0.2,temp=293,expRg=1.38,expRgErr=0.07,pH=7.5,fasta=list(fasta_Hst5),ionic=0.15,path='Hst5'+outdir)
    proteins.loc['Hst52'] = dict(eps_factor=0.2,temp=298,expRg=1.87,expRgErr=0.07,pH=7.0,fasta=list(fasta_Hst52),ionic=0.15,path='Hst52'+outdir)
    proteins.loc['aSyn140'] = dict(eps_factor=0.2,temp=293,expRg=3.55,expRgErr=0.1,pH=7.4,fasta=list(fasta_aSyn140),ionic=0.2,path='aSyn140'+outdir)
    proteins.loc['ACTR'] = dict(eps_factor=0.2,temp=278,expRg=2.63,expRgErr=0.1,pH=7.4,fasta=list(fasta_ACTR),ionic=0.2,path='ACTR'+outdir)
    proteins.loc['RNaseA'] = dict(eps_factor=0.2,temp=298,expRg=3.36,expRgErr=0.1,pH=7.5,fasta=list(fasta_RNaseA),ionic=0.15,path='RNaseA'+outdir)
    proteins.loc['PNt'] = dict(eps_factor=0.2,temp=298,expRg=5.11,expRgErr=0.2,pH=7.5,fasta=list(fasta_PNt),ionic=0.15,path='PNt'+outdir)
    proteins.loc['p15PAF'] = dict(eps_factor=0.2,temp=298,expRg=2.81,expRgErr=0.1,pH=7.0,fasta=list(fasta_p15PAF),ionic=0.15,path='p15PAF'+outdir)
    proteins.loc['CoRNID'] = dict(eps_factor=0.2,temp=293,expRg=4.7,expRgErr=0.2,pH=7.5,fasta=list(fasta_CoRNID),ionic=0.2,path='CoRNID'+outdir)
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


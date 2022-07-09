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
    return chi2 / len(prot.labels)

def optTauC(prot):
    obs = 1 if prot.obs=='ratio' else 2
    chi2list = []
    tau_c = np.arange(2,10.05,1)
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
        chi2list.append(chi2 / len(prot.labels))

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

def reweightRg(df,name,prot):
    #rg = np.sqrt( np.dot(np.power(prot.rgarray,2), prot.weights) )
    rg = np.dot(prot.rgarray, prot.weights)
    chi2_rg = np.power((prot.expRg-rg)/prot.expRgErr,2)
    #chi2_rg = np.power((prot.expRg-rg)/(prot.expRg*0.03),2)
    return rg, chi2_rg

def calcRg(df,name,prot):
    t = md.load_dcd(prot.path+"/{:s}.dcd".format(name),prot.path+"/{:s}.pdb".format(name))
    residues = [res.name for res in t.top.atoms]
    masses = df.loc[residues,'MW'].values
    masses[0] += 2
    masses[-1] += 16
    # calculate the center of mass
    cm = np.sum(t.xyz*masses[np.newaxis,:,np.newaxis],axis=1)/masses.sum()
    # calculate residue-cm distances
    si = np.linalg.norm(t.xyz - cm[:,np.newaxis,:],axis=2)
    # calculate rg
    rgarray = np.sqrt(np.sum(si**2*masses,axis=1)/masses.sum())
    #rg = np.sqrt( np.power(rgarray, 2).mean() )
    rg = rgarray.mean()
    chi2_rg = np.power((prot.expRg-rg)/prot.expRgErr,2)
    #chi2_rg = np.power((prot.expRg-rg)/(prot.expRg*0.03),2)
    return rgarray, rg, chi2_rg

def initProteins(cycle):
    outdir = '/{:d}'.format(cycle)
    proteins = pd.DataFrame(columns=['labels','wh','tau_c','temp','obs','pH','ionic','expPREs','initPREs','eff','chi2_pre','fasta','weights','path'],dtype=object)
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
    fasta_aSyn = """MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTK
EQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDP
DNEAYEMPSEEGYQDYEPEA""".replace('\n', '')
    fasta_A2 = """GHMGRGGNFGFGDSRGGGGNFGPGPGSNFRGGSDGYGSGR
GFGDGYNGYGGGPGGGNFGGSPGYGGGRGGYGGGGPGYGNQGGGYGGGYDNYGGGNYGSG
NYNDFGNYNQQPSNYGPMKSGNFGGSRNMGGPYGGGNYGPGGSGGSGGYGGRSRY""".replace('\n', '')

    proteins.loc['A2'] = dict(labels=[99, 143], tau_c=10.0,
                               wh=850,temp=298,obs='ratio',pH=5.5,fasta=list(fasta_A2),ionic=0.005,weights=False,path='A2'+outdir)
    proteins.loc['aSyn'] = dict(labels=[24, 42, 62, 87, 103], tau_c=1.0,
                               wh=700,temp=283,obs='ratio',pH=7.4,fasta=list(fasta_aSyn),ionic=0.2,weights=False,path='aSyn'+outdir)
    proteins.loc['OPN'] = dict(labels=[10, 33, 64, 88, 117, 130, 144, 162, 184, 203], tau_c=3.0,
                               wh=800,temp=298,obs='rate',pH=6.5,fasta=list(fasta_OPN),ionic=0.15,weights=False,path='OPN'+outdir)
    proteins.loc['FUS'] = dict(labels=[16, 86, 142], tau_c=10.0,
                               wh=850,temp=298,obs='rate',pH=5.5,fasta=list(fasta_FUS),ionic=0.15,weights=False,path='FUS'+outdir)
    proteins.loc['FUS12E'] = dict(labels=[16, 86, 142], tau_c=10.0,
                               wh=850,temp=298,obs='rate',pH=5.5,fasta=list(fasta_FUS12E),ionic=0.15,weights=False,path='FUS12E'+outdir)
    return proteins

def initProteinsRgs(cycle):
    outdir = '/{:d}'.format(cycle)
    proteins = pd.DataFrame(columns=['temp','expRg','expRgErr','Rg','rgarray','eff','chi2_rg','weights','pH','ionic','fasta','path'],dtype=object)
    fasta_DSS1 = """MSRAALPSLENLEDDDEFEDFATENWPMKDTELDTGDDTLWENNWDDEDIGDDDFSVQLQAELKKKGVAAC""".replace('\n', '')
    fasta_PTMA = """GPSDAAVDTSSEITTKDLKEKKEVVEEAENGRDAPANGNANEENGEQEADNEVDEEEEEGGEEEEEEEEG
DGEEEDGDEDEEAESATGKRAAEDDEDDDVDTKKQKTDEDD""".replace('\n', '')
    fasta_NHE6cmdd = """GPPLTTTLPACCGPIARCLTSPQAYENQEQLKDDDSDLILNDGDISLTYGDSTVNTEPATSSAPRRFMGNSSED
ALDRELAFGDHELVIRGTRLVLPMDDSEPPLNLLDNTRHGPA""".replace('\n', '')
    fasta_ANAC046 = """NAPSTTITTTKQLSRIDSLDNIDHLLDFSSLPPLIDPGFLGQPGPSFSGARQQHDLKPVLHHPTTAPVDNTYLPTQALNFPYHS
VHNSGSDFGYGAGSGNNNKGMIKLEHSLVSVSQETGLSSDVNTTATPEISSYPMMMNPAMMDGSKSACDGLDDLIFWEDLYTS""".replace('\n', '')
    fasta_GHRICD = """SKQQRIKMLILPPVPVPKIKGIDPDLLKEGKLEEVNTILAIHDSYKPEFHSDDSWVEFIELDIDEPDEKTEESDTDRLLSSDHEKSHSNL
GVKDGDSGRTSCCEPDILETDFNANDIHEGTSEVAQPQRLKGEADLLCLDQKNQNNSPYHDACPATQQPSVIQAEKNKPQPLPTEGAESTHQAAH
IQLSNPSSLSNIDFYAQVSDITPAGSVVLSPGQKNKAGMSQCDMHPEMVSLCQENFLMDNAYFCEADAKKCIPVAPHIKVESHIQPSLNQEDIYI
TTESLTTAAGRPGTGEHVPGSEMPVPDYTSIHIVQSPQGLILNATALPLPDKEFLSSCGYVSTDQLNKIMP""".replace('\n', '')
    fasta_Ash1 = """SASSSPSPSTPTKSGKMRSRSSSPVRPKAYTPSPRSPNYHRFALDSPPQSPRRSSNSSITKKGSRRSSGSSPTRHTTRVCV"""
    fasta_CTD2 = """FAGSGSNIYSPGNAYSPSSSNYSPNSPSYSPTSPSYSPSSPSYSPTSPCYSPTSPSYSPTSPNYTPVTPSYSPTSPNYSASPQ"""
    fasta_Hst52 = """DSHAKRHHGYKRKFHEKHHSHRGYDSHAKRHHGYKRKFHEKHHSHRGY""".replace('\n', '') # DOI: 10.1021/acs.jpcb.0c09635
    fasta_Hst5 = """DSHAKRHHGYKRKFHEKHHSHRGY""".replace('\n', '')
    fasta_aSyn140 = """MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTK
EQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA""".replace('\n', '')
    fasta_PNt = """DWNNQSIVKTGERQHGIHIQGSDPGGVRTASGTTIKVSGRQAQGILLENPAAELQFRNGSVTSS
GQLSDDGIRRFLGTVTVKAGKLVADHATLANVGDTWDDDGIALYVAGEQAQASIADSTLQGAGG
VQIERGANVTVQRSAIVDGGLHIGALQSLQPEDLPPSRVVLRDTNVTAVPASGAPAAVSVLGAS
ELTLDGGHITGGRAAGVAAMQGAVVHLQRATIRRGEALAGGAVPGGAVPGGAVPGGFGPGGFGP
VLDGWYGVDVSGSSVELAQSIVEAPELGAAIRVGRGARVTVPGGSLSAPHGNVIETGGARRFAP
QAAPLSITLQAGAH""".replace('\n','')
    fasta_PNtS1 = """DWNNQSIVKTGERQHGIHIQGSDPGGVRTASGTTIKVSGRQAQGILLENPAAELQFRNGSVTSS
GQKSDDGIRRFLGTVTVLAGKLVADHATLANVGDTWDDDGIALYVAGEQAQASIADSTLQGAGG
VQIERGANVTVQRSAIVLGGLHIGALQSLQPEDDPPSRVVLRDTNVTAVPASGAPAAVSVLGAS
LLTLDGGHITGGRAAGVAAMQGAVVHEQRATIRRGEALAGGAVPGGAVPGGAVPGGFGPGGFGP
VLDGWYGVDVSGSSVELAQSIVEAPELGAAIRVGRGARVTVPGGSLSAPHGNVIETGGARRFAP
QAAPLSITLQAGAH""".replace('\n','')
    fasta_PNtS4 = """DWNNQSIVKTGERQHGIHIQGSDPGGVRTASGTTIKVSGRQAQGILLENPAAELQFRNGSVTSS
GQLSFVGITRDLGRDTVKAGKLVADHATLANVGDTWDDDGIALYVAGEQAQASIADSTLQGAGG
VQIERGADVRVQREAIVDGGLHNGALQSLQPSILPPSTVVLRDTNVTAVPASGAPAAVLVSGAS
GLRLDGGHIHEGRAAGVAAMQGAVVTLQTATIRRGEALAGGAVPGGAVPGGAVPGGFGPGGFGP
VLDGWYGVDVSGSSVELAQSIVEAPELGAAIRVGRGARVTVPGGSLSAPHGNVIETGGARRFAP
QAAPLSITLQAGAH""".replace('\n','')
    fasta_PNtS5 = """DWNNQSIVKTGERQHGIHIQGSDPGGVRTASGTTIKVSGRQAQGILLENPAAELQFRNGSVTSS
GQLSDDGIEDFLGTVTVDAGELVADHATLANVGDTWDDDGIALYVAGEQAQASIADSTLQGAGG
VQIEDGANVTVQESAIVDGGLHIGALQSLQPRRLPPSRVVLRKTNVTAVPASGAPAAVSVLGAS
KLTLRGGHITGGRAAGVAAMQGAVVHLQRATIRRGRALAGGAVPGGAVPGGAVPGGFGPGGFGP
VLDGWYGVDVSGSSVELAQSIVEAPELGAAIRVGRGARVTVPGGSLSAPHGNVIETGGARRFAP
QAAPLSITLQAGAH""".replace('\n','')
    fasta_PNtS6 = """DWNNQSIVKTGERQHGIHIQGSDPGGVRTASGTTIKVSGRQAQGILLENPAAELQFRNGSVTSS
GQLSDRGIDRFLGTVTVEAGKLVADHATLANVGDTWDKDGIALYVAGRQAQASIADSTLQGAGG
VQIREGANVTVQRSAIVDGGLHIGALQSLQPERLPPSDVVLRDTNVTAVPASGAPAAVSVLGAS
RLTLDGGHITGGDAAGVAAMQGAVVHLQRATIERGEALAGGAVPGGAVPGGAVPGGFGPGGFGP
VLDGWYGVDVSGSSVELAQSIVEAPELGAAIRVGRGARVTVPGGSLSAPHGNVIETGGARRFAP
QAAPLSITLQAGAH""".replace('\n','')
    fasta_ACTR = """GTQNRPLLRNSLDDLVGPPSNLEGQSDERALLDQLHTLLSNTDATGLEEIDRALGIPELVNQGQALEPKQD""".replace('\n', '') # DOI: 10.1073/pnas.1322611111
    fasta_RNaseA = """KETAAAKFERQHMDSSTSAASSSNYCNQMMKSRNLTKDRCKPVNTFVHESLADVQAVCSQKNVACKNGQTNCYQSYSTMSI
TDCRETGSSKYPNCAYKTTQANKHIIVACEGNPYVPVHFDASV""".replace('\n', '')
    fasta_p15PAF = """MVRTKADSVPGTYRKVVAARAPRKVLGSSTSATNSTSVSSRKAENKYAGGNPVCVRPTPKWQKGIGEFFR
LSPKDSEKENQIPEEAGSSGLGKAKRKACPLQPDHTNDEKE""".replace('\n', '') # DOI: 10.1016/j.bpj.2013.12.046
    fasta_CoRNID = """GPHMQVPRTHRLITLADHICQIITQDFARNQVPSQASTSTFQTSPSALSSTPVRTKTSSRYS
PESQSQTVLHPRPGPRVSPENLVDKSRGSRPGKSPERSHIPSEPYEPISPPQGPAVHEKQDSMLLLSQRGVDPAEQRSDSRSP
GSISYLPSFFTKLESTSPMVKSKKQEIFRKLNSSGGGDSDMAAAQPGTEIFNLPAVTTSGAVSSRSHSFADPASNLGLEDIIR
KALMGSFDDKVEDHGVVMSHPVGIMPGSASTSVVTSSEARRDE""".replace('\n', '') # SASDF34
    fasta_Sic1 = """GSMTPSTPPRSRGTRYLAQPSGNTSSSALMQGQKTPQKPS
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
    fasta_tau35 = """EPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAP
VPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGS
VQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTD
HGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL""".replace('\n', '')
    fasta_CAHSD = """MSGRNVESHMERNEKVVVNNSGHADVKKQQQQVEHTEFTHTEVKAPLIHPAPPIISTGAAGLA
EEIVGQGFTASAARISGGTAEVHLQPSAAMTEEARRDQERYRQEQESIAKQQEREMEKKTEAYRKT
AEAEAEKIRKELEKQHARDVEFRKDLIESTIDRQKREVDLEAKMAKRELDREGQLAKEALERSRLA
TNVEVNFDSAAGHTVSGGTTVSTSDKMEIKRN""".replace('\n','')
    fasta_p532070 = """GPGSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAALEHHHHHH""" # DOI: 10.1038/s41467-021-21258-5
    fasta_ht2N3R = """MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPG
SETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAG
HVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPK
TPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAK
SRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIVYKPVDLSKVTSKCGSLGNIHHK
PGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAE
IVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL""".replace('\n','')
    fasta_ht2N4R = """MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPG
SETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAG
HVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPK
TPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAK
SRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHV
PGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNI
THVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMV
DSPQLATLADEVSASLAKQGL""".replace('\n', '')
    fasta_p27Cv14 = """GSHMKGACKSSSPPSNDQGRPGDPKQVIDKTEVERTQDTSNIQETQSANNSGPDKPS
RCDLAVSGVAAAALPAPGHANSTARDLTRDEEAGSVEQTPKKPGLRRRQT""".replace('\n', '')
    fasta_p27Cv15 = """GSHMKGACIVANSPPDDVKSKEDVPQTDPRLTGGDRDNARASRTGNDPAGASTQSAE
VACSNPILSTPDAQEKQAGTSNSKERPHEQLSAGSVEQTPKKPGLRRRQT""".replace('\n', '')
    fasta_p27Cv31 = """GSHMKGACKVPAQESQDVSGSRPAAPLIGAPANSEDTHLVDPKTDPSDSQTGLAEQC
AGIRKRPATDDSSTQNKRANRTEENVSDGSPNAGSVEQTPKKPGLRRRQT""".replace('\n', '')
    fasta_p27Cv44 = """GSHMKGACRKPANAEADSSSCQNVPRGKSKQAPETPTGSPLGDATLNQVKPRRPSSA
STNIGQLEDADEDDAEDHVGSAVTSQTIPNDRAGSVEQTPKKPGLRRRQT""".replace('\n', '')
    fasta_p27Cv56 = """GSHMKGACGSSVLGTGNPRNQAHVSDTSLEEDDDEQDDSTPDEVSQACTIVASALDI
NAATPRSPKASPKRKRKRQSTAPAQGNEPPGNAGSVEQTPKKPGLRRRQT""".replace('\n', '')
    fasta_p27Cv78 = """GSHMKGACALPSGVVPAEDDDDDEEEEDDQDPAQPQAVQGAAPSSGTNNSQPILPSIA
VNSTTGPNSTAGKKKRKRRRTRHSNCATLSSAGSVEQTPKKPGLRRRQT""".replace('\n', '')
    fasta_A1S = """GSMASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGYGGSGDGYNGFGNDGSNFGGGGNYNNQ
SSNFGPMKGGNFGGRSSGPYGGGGQYFAKPRNQGGYGGSSSSSSYGSGRRF""".replace('\n', '')

    proteins.loc['A1S50'] = dict(temp=293,expRg=2.645,expRgErr=0.02,pH=7.5,fasta=list(fasta_A1S),ionic=0.05,path='A1S50'+outdir)
    proteins.loc['A1S150'] = dict(temp=293,expRg=2.65,expRgErr=0.02,pH=7.5,fasta=list(fasta_A1S),ionic=0.15,path='A1S150'+outdir)
    proteins.loc['A1S300'] = dict(temp=293,expRg=2.62,expRgErr=0.02,pH=7.5,fasta=list(fasta_A1S),ionic=0.3,path='A1S300'+outdir)
    proteins.loc['A1S500'] = dict(temp=293,expRg=2.528,expRgErr=0.02,pH=7.5,fasta=list(fasta_A1S),ionic=0.5,path='A1S500'+outdir)
    proteins.loc['p27Cv14'] = dict(temp=293,expRg=2.936,expRgErr=0.13,pH=7.2,fasta=list(fasta_p27Cv14),ionic=0.095,path='p27Cv14'+outdir)
    proteins.loc['p27Cv15'] = dict(temp=293,expRg=2.915,expRgErr=0.1,pH=7.2,fasta=list(fasta_p27Cv15),ionic=0.095,path='p27Cv15'+outdir)
    proteins.loc['p27Cv31'] = dict(temp=293,expRg=2.81,expRgErr=0.18,pH=7.2,fasta=list(fasta_p27Cv31),ionic=0.095,path='p27Cv31'+outdir)
    proteins.loc['p27Cv44'] = dict(temp=293,expRg=2.492,expRgErr=0.13,pH=7.2,fasta=list(fasta_p27Cv44),ionic=0.095,path='p27Cv44'+outdir)
    proteins.loc['p27Cv56'] = dict(temp=293,expRg=2.328,expRgErr=0.1,pH=7.2,fasta=list(fasta_p27Cv56),ionic=0.095,path='p27Cv56'+outdir)
    proteins.loc['p27Cv78'] = dict(temp=293,expRg=2.211,expRgErr=0.03,pH=7.2,fasta=list(fasta_p27Cv78),ionic=0.095,path='p27Cv78'+outdir)
    proteins.loc['tau35'] = dict(temp=298,expRg=4.64,expRgErr=0.1,pH=7.4,fasta=list(fasta_tau35),ionic=0.15,path='tau35'+outdir)
    proteins.loc['ht2N3R'] = dict(temp=298,expRg=6.3,expRgErr=0.3,pH=7.4,fasta=list(fasta_ht2N3R),ionic=0.15,path='ht2N3R'+outdir)
    proteins.loc['ht2N4R'] = dict(temp=298,expRg=6.7,expRgErr=0.3,pH=7.4,fasta=list(fasta_ht2N4R),ionic=0.15,path='ht2N4R'+outdir)
    proteins.loc['CAHSD'] = dict(temp=293,expRg=4.84,expRgErr=0.2,pH=7.0,fasta=list(fasta_CAHSD),ionic=0.07,path='CAHSD'+outdir)
    proteins.loc['DSS1'] = dict(temp=288,expRg=2.5,expRgErr=0.1,pH=7.4,fasta=list(fasta_DSS1),ionic=0.17,path='DSS1'+outdir)
    proteins.loc['PTMA'] = dict(temp=288,expRg=3.7,expRgErr=0.2,pH=7.4,fasta=list(fasta_PTMA),ionic=0.16,path='PTMA'+outdir)
    proteins.loc['NHE6cmdd'] = dict(temp=288,expRg=3.2,expRgErr=0.2,pH=7.4,fasta=list(fasta_NHE6cmdd),ionic=0.17,path='NHE6cmdd'+outdir)
    proteins.loc['ANAC046'] = dict(temp=298,expRg=3.6,expRgErr=0.3,pH=7.0,fasta=list(fasta_ANAC046),ionic=0.14,path='ANAC046'+outdir)
    proteins.loc['GHRICD'] = dict(temp=298,expRg=6.0,expRgErr=0.5,pH=7.3,fasta=list(fasta_GHRICD),ionic=0.35,path='GHRICD'+outdir)
    proteins.loc['p532070'] = dict(eps_factor=0.2,temp=277,expRg=2.39,expRgErr=0.05,pH=7,fasta=list(fasta_p532070),ionic=0.1,path='p532070'+outdir)
    proteins.loc['Ash1'] = dict(temp=293,expRg=2.9,expRgErr=0.05,pH=7.5,fasta=list(fasta_Ash1),ionic=0.15,path='Ash1'+outdir)
    proteins.loc['CTD2'] = dict(temp=293,expRg=2.614,expRgErr=0.05,pH=7.5,fasta=list(fasta_CTD2),ionic=0.12,path='CTD2'+outdir)
    proteins.loc['ColNT'] = dict(temp=277,expRg=2.83,expRgErr=0.1,pH=7.6,fasta=list(fasta_ColNT),ionic=0.4,path='ColNT'+outdir)
    proteins.loc['hNL3cyt'] = dict(temp=293,expRg=3.15,expRgErr=0.2,pH=8.5,fasta=list(fasta_hNL3cyt),ionic=0.3,path='hNL3cyt'+outdir)
    proteins.loc['SH4UD'] = dict(temp=293,expRg=2.71,expRgErr=0.1,pH=8.0,fasta=list(fasta_SH4UD),ionic=0.2,path='SH4UD'+outdir)
    proteins.loc['Sic1'] = dict(temp=293,expRg=3.0,expRgErr=0.4,pH=7.5,fasta=list(fasta_Sic1),ionic=0.2,path='Sic1'+outdir)
    proteins.loc['FhuA'] = dict(temp=298,expRg=3.34,expRgErr=0.1,pH=7.5,fasta=list(fasta_FhuA),ionic=0.15,path='FhuA'+outdir)
    proteins.loc['K10'] = dict(temp=288,expRg=4.0,expRgErr=0.1,pH=7.4,fasta=list(fasta_K10),ionic=0.15,path='K10'+outdir)
    proteins.loc['K27'] = dict(temp=288,expRg=3.7,expRgErr=0.2,pH=7.4,fasta=list(fasta_K27),ionic=0.15,path='K27'+outdir)
    proteins.loc['K25'] = dict(temp=288,expRg=4.1,expRgErr=0.2,pH=7.4,fasta=list(fasta_K25),ionic=0.15,path='K25'+outdir)
    proteins.loc['K32'] = dict(temp=288,expRg=4.2,expRgErr=0.3,pH=7.4,fasta=list(fasta_K32),ionic=0.15,path='K32'+outdir)
    proteins.loc['K23'] = dict(temp=288,expRg=4.9,expRgErr=0.2,pH=7.4,fasta=list(fasta_K23),ionic=0.15,path='K23'+outdir)
    proteins.loc['K44'] = dict(temp=288,expRg=5.2,expRgErr=0.2,pH=7.4,fasta=list(fasta_K44),ionic=0.15,path='K44'+outdir)
    proteins.loc['A1'] = dict(temp=298,expRg=2.76,expRgErr=0.02,pH=7.0,fasta=list(fasta_A1),ionic=0.15,path='A1'+outdir)
    proteins.loc['M12FP12Y'] = dict(temp=298,expRg=2.60,expRgErr=0.02,pH=7.0,fasta=list(fasta_M12FP12Y),ionic=0.15,path='M12FP12Y'+outdir)
    proteins.loc['P7FM7Y'] = dict(temp=298,expRg=2.72,expRgErr=0.01,pH=7.0,fasta=list(fasta_P7FM7Y),ionic=0.15,path='P7FM7Y'+outdir)
    proteins.loc['M9FP6Y'] = dict(temp=298,expRg=2.66,expRgErr=0.01,pH=7.0,fasta=list(fasta_M9FP6Y),ionic=0.15,path='M9FP6Y'+outdir)
    proteins.loc['M8FP4Y'] = dict(temp=298,expRg=2.71,expRgErr=0.01,pH=7.0,fasta=list(fasta_M8FP4Y),ionic=0.15,path='M8FP4Y'+outdir)
    proteins.loc['M9FP3Y'] = dict(temp=298,expRg=2.68,expRgErr=0.01,pH=7.0,fasta=list(fasta_M9FP3Y),ionic=0.15,path='M9FP3Y'+outdir)
    proteins.loc['M10R'] = dict(temp=298,expRg=2.67,expRgErr=0.01,pH=7.0,fasta=list(fasta_M10R),ionic=0.15,path='M10R'+outdir)
    proteins.loc['M6R'] = dict(temp=298,expRg=2.57,expRgErr=0.01,pH=7.0,fasta=list(fasta_M6R),ionic=0.15,path='M6R'+outdir)
    proteins.loc['P2R'] = dict(temp=298,expRg=2.62,expRgErr=0.02,pH=7.0,fasta=list(fasta_P2R),ionic=0.15,path='P2R'+outdir)
    proteins.loc['P7R'] = dict(temp=298,expRg=2.71,expRgErr=0.01,pH=7.0,fasta=list(fasta_P7R),ionic=0.15,path='P7R'+outdir)
    proteins.loc['M3RP3K'] = dict(temp=298,expRg=2.63,expRgErr=0.02,pH=7.0,fasta=list(fasta_M3RP3K),ionic=0.15,path='M3RP3K'+outdir)
    proteins.loc['M6RP6K'] = dict(temp=298,expRg=2.79,expRgErr=0.01,pH=7.0,fasta=list(fasta_M6RP6K),ionic=0.15,path='M6RP6K'+outdir)
    proteins.loc['M10RP10K'] = dict(temp=298,expRg=2.85,expRgErr=0.01,pH=7.0,fasta=list(fasta_M10RP10K),ionic=0.15,path='M10RP10K'+outdir)
    proteins.loc['M4D'] = dict(temp=298,expRg=2.64,expRgErr=0.01,pH=7.0,fasta=list(fasta_M4D),ionic=0.15,path='M4D'+outdir)
    proteins.loc['P4D'] = dict(temp=298,expRg=2.72,expRgErr=0.03,pH=7.0,fasta=list(fasta_P4D),ionic=0.15,path='P4D'+outdir)
    proteins.loc['P8D'] = dict(temp=298,expRg=2.69,expRgErr=0.01,pH=7.0,fasta=list(fasta_P8D),ionic=0.15,path='P8D'+outdir)
    proteins.loc['P12D'] = dict(temp=298,expRg=2.80,expRgErr=0.01,pH=7.0,fasta=list(fasta_P12D),ionic=0.15,path='P12D'+outdir)
    proteins.loc['P12E'] = dict(temp=298,expRg=2.85,expRgErr=0.01,pH=7.0,fasta=list(fasta_P12E),ionic=0.15,path='P12E'+outdir)
    proteins.loc['P7KP12D'] = dict(temp=298,expRg=2.92,expRgErr=0.01,pH=7.0,fasta=list(fasta_P7KP12D),ionic=0.15,path='P7KP12D'+outdir)
    proteins.loc['P7KP12Db'] = dict(temp=298,expRg=2.56,expRgErr=0.01,pH=7.0,fasta=list(fasta_P7KP12Db),ionic=0.15,path='P7KP12Db'+outdir)
    proteins.loc['M12FP12YM10R'] = dict(temp=298,expRg=2.61,expRgErr=0.02,pH=7.0,fasta=list(fasta_M12FP12YM10R),ionic=0.15,path='M12FP12YM10R'+outdir)
    proteins.loc['M10FP7RP12D'] = dict(temp=298,expRg=2.86,expRgErr=0.01,pH=7.0,fasta=list(fasta_M10FP7RP12D),ionic=0.15,path='M10FP7RP12D'+outdir)
    proteins.loc['Hst5'] = dict(temp=293,expRg=1.38,expRgErr=0.05,pH=7.5,fasta=list(fasta_Hst5),ionic=0.15,path='Hst5'+outdir)
    proteins.loc['Hst52'] = dict(temp=298,expRg=1.87,expRgErr=0.05,pH=7.0,fasta=list(fasta_Hst52),ionic=0.15,path='Hst52'+outdir)
    proteins.loc['aSyn140'] = dict(temp=293,expRg=3.55,expRgErr=0.1,pH=7.4,fasta=list(fasta_aSyn140),ionic=0.2,path='aSyn140'+outdir)
    proteins.loc['ACTR'] = dict(temp=278,expRg=2.63,expRgErr=0.1,pH=7.4,fasta=list(fasta_ACTR),ionic=0.2,path='ACTR'+outdir)
    proteins.loc['RNaseA'] = dict(temp=298,expRg=3.36,expRgErr=0.1,pH=7.5,fasta=list(fasta_RNaseA),ionic=0.15,path='RNaseA'+outdir)
    proteins.loc['p15PAF'] = dict(temp=298,expRg=2.81,expRgErr=0.1,pH=7.0,fasta=list(fasta_p15PAF),ionic=0.15,path='p15PAF'+outdir)
    proteins.loc['CoRNID'] = dict(temp=293,expRg=4.7,expRgErr=0.2,pH=7.5,fasta=list(fasta_CoRNID),ionic=0.2,path='CoRNID'+outdir)
    proteins.loc['PNt'] = dict(temp=298,expRg=5.11,expRgErr=0.1,pH=7.5,fasta=list(fasta_PNt),ionic=0.15,path='PNt'+outdir)
    proteins.loc['PNtS1'] = dict(temp=298,expRg=4.92,expRgErr=0.1,pH=7.5,fasta=list(fasta_PNtS1),ionic=0.15,path='PNtS1'+outdir)
    proteins.loc['PNtS4'] = dict(temp=298,expRg=5.34,expRgErr=0.1,pH=7.5,fasta=list(fasta_PNtS4),ionic=0.15,path='PNtS4'+outdir)
    proteins.loc['PNtS5'] = dict(temp=298,expRg=4.87,expRgErr=0.1,pH=7.5,fasta=list(fasta_PNtS5),ionic=0.15,path='PNtS5'+outdir)
    proteins.loc['PNtS6'] = dict(temp=298,expRg=5.26,expRgErr=0.1,pH=7.5,fasta=list(fasta_PNtS6),ionic=0.15,path='PNtS6'+outdir)
    return proteins

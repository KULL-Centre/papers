import pandas as pd
import numpy as np
import mdtraj as md
import itertools

def initProteins():
    proteins = pd.DataFrame(index=['aSyn','FUS','A2'], columns=['labels','wh','temp','obs','pH','expPREs','ionic','fasta','path'])
    fasta_aSyn = """MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTK
EQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDP
DNEAYEMPSEEGYQDYEPEA""".replace('\n', '')
    fasta_FUS = """MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQ
SQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSS
SYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNS""".replace('\n', '')
    fasta_A2 = """GHMGRGGNFGFGDSRGGGGNFGPGPGSNFRGGSDGYGSGR
GFGDGYNGYGGGPGGGNFGGSPGYGGGRGGYGGGGPGYGNQGGGYGGGYDNYGGGNYGSG
NYNDFGNYNQQPSNYGPMKSGNFGGSRNMGGPYGGGNYGPGGSGGSGGYGGRSRY""".replace('\n', '')

    proteins.loc['A2'] = dict(labels=[99, 143],
                               wh=850,temp=298,obs='ratio',pH=5.5,fasta=list(fasta_A2),ionic=0.005,path='A2')
    proteins.loc['FUS'] = dict(labels=[16, 86, 142],
                               wh=850,temp=298,obs='rate',pH=5.5,fasta=list(fasta_FUS),ionic=0.15,path='FUS')
    proteins.loc['aSyn'] = dict(labels=[24, 42, 62, 87, 103],
                               wh=700,temp=283,obs='ratio',pH=7.4,fasta=list(fasta_aSyn),ionic=0.15,path='aSyn')
    return proteins


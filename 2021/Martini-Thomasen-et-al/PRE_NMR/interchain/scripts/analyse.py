import pandas as pd
import numpy as np
import mdtraj as md
import itertools

def initProteins():
    proteins = pd.DataFrame(index=['FUS'], columns=['labels','eps_factor','wh','L','temp','obs','pH','ionic','expPREs','fasta','path'])
    fasta_FUS = """MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQ
SQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSS
SYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNS""".replace('\n', '')
    proteins.loc['FUS'] = dict(labels=[16, 86, 142],eps_factor=0.2,L=40.5,wh=850,temp=298,obs='rate',pH=5.5,fasta=list(fasta_FUS),ionic=0.15,path='FUS/') 
    return proteins


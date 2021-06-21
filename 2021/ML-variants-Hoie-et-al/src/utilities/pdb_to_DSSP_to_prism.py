#!/lindorffgrp-isilon/maghoi/_conda_env/envs/DSSP/bin/python

#Tuple Index	Value
#0	DSSP index
#1	Amino acid
#2	Secondary structure
#3	Relative ASA
#4	Phi
#5	Psi
#6	NH-->O_1_relidx
#7	NH-->O_1_energy
#8	O-->NH_1_relidx
#9	O-->NH_1_energy
#10	NH-->O_2_relidx
#11	NH-->O_2_energy
#12	O-->NH_2_relidx
#13	O-->NH_2_energy
#https://biopython.org/DIST/docs/api/Bio.PDB.DSSP%27-module.html

import numpy as np
import pandas as pd
import Bio, sys, os, re, glob
from argparse import ArgumentParser
import prism_analysis as prism
prismfolder = prism.prism_folder()
import PrismData as PrismData
prismparser = PrismData.PrismParser()

def is_valid_file(parser, arg):
    if not os.path.exists(arg): parser.error("Input file %s does not exist" % arg)
    else: return open(arg, 'r')  # return an open file handle

parser = ArgumentParser(description="Convert Rosetta output ddg file to variants with score in .prism format")
parser.add_argument("-i", dest="INPUT_FILE", required=True,
                    help="input file", metavar="FILE",
                    type=lambda x: is_valid_file(parser, x))
parser.add_argument("-n", dest="PROTEIN_NAME", default="protein", required=False,
                    help="protein name")
parser.add_argument("-u", dest="UNIPROT", default="UNIPROT", required=False,
                    help="uniprot id")
parser.add_argument("-v", dest="VERBOSE", default=0, required=False,
                    help="Print error messages")
args = parser.parse_args()

INPUT_FILE = str(args.INPUT_FILE.name)
PROTEIN_NAME = str(args.PROTEIN_NAME)
UNIPROT = str(args.UNIPROT)
verbose = int(args.VERBOSE)

print("Extracting PSSM ", INPUT_FILE, "to PSSM")

# Definitions
aa_order_alphabetical = pd.Series(["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",
           "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])

# Custom pdb_dict
pdb_dict = {
    "2GRN":["UBE2I", "P63279"],
    "1WYW":["SUMO1", "P63165"],
    "3S4Y":["TPK1", "Q9H3S4"],
    "4DJC":["CALM1", "P0DP23"],
    "1JM7":["BRCA1", "P38398"],
    "6RZ3":["P53", "P04637"],
    "4Q21":["HRAS", "P01112"],
    "2OJG":["MAPK1", "P28482"],
    "1D5R":["PTEN", "P60484"],
    "2BZG":["TPMT", "P51580"],
    "3KJ6":["ADRB2", "P07550"],
    "1A4H":["HSP82", "P02829"],
    "3OLM":["UBI4", "P0CG63"],
    "6RSK":["PAB1", "P04147"],
    "3COQ":["GAL4", "P04386"],
    "1AH9":["IF-1", "P69222"],
    "6BVC":["GmR", "Q53396"],
    "1LI9":["bla", "P62593"],
    "1VUB":["ccdB", "P62554"],
    "3UBT":["haeIIIM", "P20589"],
    "5FW1":["cas9", "Q99ZW2"],
    "4UN3":["env", "P03377"],
    "6MYA":["HAh1n1", "A0A2Z5U3Z0"],
    "2YP7":["HAh3n2", "A0A097PF60"],
    "6QNW": ["PAh1n1", "P15659"],
}


### Load DSSP
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
pdbparser = PDBParser()

pdb_name = os.path.splitext(os.path.basename(INPUT_FILE))[0]
pdb_filepath = INPUT_FILE

structure = pdbparser.get_structure(pdb_name, pdb_filepath)
model = structure[0]
dssp = DSSP(model, pdb_filepath)

# Extract seq fromp pdb
# Generate empty vector to fill in seq
last_pos = dssp[dssp.keys()[-1]][0]
seq = np.zeros(last_pos, dtype = str)
seq[:] = "X"

for i, key in enumerate(list(dssp.keys())):
    d_all = dssp[key]
    d_i, d_aa, d_ss, d_rsa = d_all[0:4]
    seq[d_i-1] = d_aa
seq = "".join(seq)

# Create prism_df to fill
variant_df = prism.generate_empty_prism_df(seq)
variant_df["ss"], variant_df["rsa"] = np.nan, np.nan
dssp_df = variant_df[["variant", "ss", "rsa"]]

# Fill with DSSP values
for i, key in enumerate(list(dssp.keys())):    
    d_all = dssp[key]
    d_i, d_aa, d_ss, d_rsa = d_all[0:4]
    
    for mut in aa_order_alphabetical:
        v = str(d_aa)+str(d_i)+str(mut)
        dssp_df.loc[v] = (v, d_ss, d_rsa)
        
dssp_df = dssp_df.dropna() 

# Create prismdata container, fill with header and dataframe
#try: protein_name, uniprot = pdb_dict[pdb_name]
try: protein_name, uniprot = PROTEIN_NAME, UNIPROT
except:
    print("Cannot find", pdb_name, "in pdb_dict. Using default names")
    protein_name, uniprot = "protein", "uniprot"
basename = protein_name + "_" + uniprot + "_" + pdb_name
    
# Load "empty" prismdata to replace
prism_data = prismfolder.generate_empty_prismdata(prism_type="prism_rosetta")

# Generate new metadata dict
metadata, outfile = prismfolder.fill_prism_header(protein_name, seq, uniprot, basename, prism_type = "dssp")

# update dict and dataframe
metadata["protein"]["pdb"] = pdb_name
metadata["dssp"] = metadata.pop("gemme")
metadata["columns"] = {"ss": "Secondary structure extracted from pdb using DSSP",
              "rsa": "Relative solvent accessbility extracted from pdb using DDSP"}
dataframe = dssp_df

# more updates
from datetime import date
today = str(date.today())
metadata["date"] = today
metadata["dssp"]["date"] = today

# Assign and save
prism_data.metadata, prism_data.dataframe = metadata, dataframe
print("Writing prismfile to", outfile)
prismparser.write(outfile, prism_data)

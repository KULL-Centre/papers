#!/lindorffgrp-isilon/maghoi/_conda_env/bin/python
import pandas as pd
import numpy as np
import re

from argparse import ArgumentParser
import os
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("Input file %s does not exist" % arg)
    else:
        return open(arg, 'r')  # return an open file handle

parser = ArgumentParser(description="Convert MAVEdb file to .prism format")
parser.add_argument("-i", dest="INPUT_FILE", required=True,
                    help="input file", metavar="FILE",
                    type=lambda x: is_valid_file(parser, x))
parser.add_argument("-o", dest="OUTPUT_FILE",
                    help="output file", metavar="FILE",
                    type=lambda x: is_valid_file(parser, x))
parser.add_argument("-c", dest="COLUMNS", default = "['score']",
                    help="Columns, string python format")
args = parser.parse_args()

INFILE = str(args.INPUT_FILE.name)
try: OUTFILE = str(args.OUTPUT_FILE.name)
except: OUTFILE = False

try: COLUMNS = eval(args.COLUMNS)
except: 
    print("Unable to extract columns. Using python list format with strings?")
    import sys
    sys.exit()

# HVGS encoding
# http://www.hgvs.org/mutnomen/references.html
aa_letter_to_threeletter = {"A": "Ala",
"V": "Val",
"I": "Ile",
"L": "Leu",
"M": "Met",
"F": "Phe",
"Y": "Tyr",
"W": "Trp",
"S": "Ser",
"T": "Thr",
"N": "Asn",
"Q": "Gln",
"C": "Cys",
"G": "Gly",
"P": "Pro",
"R": "Arg",
"H": "His",
"K": "Lys",
"D": "Asp",
"E": "Glu",
"B": "Asx",
"U": "Sec",
"X": "Xaa",
"Z": "Glx",
"*": "*", #Termination
"*": "Ter", # Termination
"~": "del", # Deletion
}

aa_threeletter_to_letter = {v: k for k, v in aa_letter_to_threeletter.items()}

def mave_to_variant(text, v = False):
    # Captures all 3 letter + 1-4 digits + 3 letters
    # Also if last (non-capturing) group is =
    p = r'([A-Za-z]{3}\d{1,4}(?:[A-Za-z]{3}|=))'
    n_mutant_text = re.findall(p, text)
    if v: print(text); print(n_mutant_text)
    
    n_mutations_list = []
    for single_mutation_text in n_mutant_text:
        p2 = r'([A-Za-z]{3})(\d{1,4})([A-Za-z]{3}|=)'
        single_mutation_re = re.match(p2, single_mutation_text)
        
        try: wt_re, pos_re, mut_re = single_mutation_re.groups()
        except:
            if text == "_wt":
                return("WT")
            else:
                print("Failed extract", text)
                return("FAILED")
        if v: print(wt_re, pos_re, mut_re)

        try: pos = str(int(pos_re))
        except: print("Pos failed", pos_re)

        try: wt = aa_threeletter_to_letter[wt_re]
        except: print("wt failed", wt_re)

        try: mut = aa_threeletter_to_letter[mut_re]
        except:
            if mut_re == "=" and wt:
                mut = wt
            else: 
                print("mut failed", mut_re)
                
        if v: print(wt, pos, mut)
            
        # Merge as string, append to n mutant list
        n_mutant = "".join([wt, pos, mut])
        n_mutations_list.append(n_mutant)
        if len(n_mutant) > 5 or len(n_mutant) < 3:
            print("Unusual length (?):", n_mutant)
        if "X" in n_mutant or "B" in n_mutant or "Z" in n_mutant:
            print("Non-specific residue (?):", n_mutant)
    
    variant = ":".join(n_mutations_list)

    #print(variant)
    return(variant)

# Read file
df_raw = pd.read_csv(INFILE, skiprows=4)
print("Shape:", df_raw.shape)
print("Columns:", df_raw.columns)
print("Statistics:\n", df_raw.describe())

# HGVS_pro to wt pos mut
#variants = df_raw["hgvs_pro"].apply(lambda x: "".join(mave_to_variant(x)))

# Check for custom extract columns
cols = ["hgvs_pro"] + COLUMNS
keep_cols = ["variant"] + COLUMNS
print("Extracting", cols)

df_proc = df_raw[cols]
variants = df_raw["hgvs_pro"].apply(lambda x: "".join(mave_to_variant(x)))
df_proc.insert(0, "variant", variants)
df_proc = df_proc[keep_cols]

# Put WT row first
wt_row_exclude = df_proc[df_proc["variant"] == "WT"].index
wt_row = df_proc.loc[wt_row_exclude]
df_proc = df_proc.drop(wt_row_exclude, axis = 0)
df_proc = pd.concat([wt_row, df_proc])

# Remove missing
missing = df_proc["variant"] == "FAILED"
df_proc = df_proc[~missing]

# Keep only specific columns
df_proc = df_proc[keep_cols]

output_filename = INFILE + ".txt"
print("\Saving file to:", output_filename)
df_proc.to_csv(output_filename, sep = "\t", index = False, na_rep=np.NaN)
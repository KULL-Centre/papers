import os, sys, glob, re
import pandas as pd
import numpy as np
import PrismData as PrismData

verbose = True
prismparser = PrismData.PrismParser()
outdir = "data/preprocessed/"

def get_filepath_dict(folder, verbose = 1):
    filepaths = glob.glob(folder + "/prism_*.txt")

    from collections import defaultdict
    protein_name_filepath_dict = defaultdict(list)

    for file in filepaths:
        file_name = os.path.basename(file)
        z = re.findall("([a-zA-Z0-9-]{3,10})[_\.]", file_name)
        try:
            gene = z[3]
            protein_name_filepath_dict[gene].append(file)
            #print(gene, uniprot)
        except:
            if verbose > 0: print("Unable to extract:", file_name)
    return(protein_name_filepath_dict)

# Load data from prism_folder
mave_path_dict = get_filepath_dict(folder = "data/raw/maves")
gemme_path_dict = get_filepath_dict(folder = "data/raw/gemme")
rosetta_path_dict = get_filepath_dict(folder = "data/raw/rosetta")
dssp_path_dict = get_filepath_dict(folder = "data/raw/dssp")

# Assign to path dict list
path_dict_list = [mave_path_dict, gemme_path_dict, rosetta_path_dict, dssp_path_dict]
file_count = len(path_dict_list)

# Generate merged prism files in folder (must exist)
os.makedirs(outdir, exist_ok = True)

# Remove multivariants
def remove_multivariants_prismdata(prismdata):
    r = r'^[ACDEFGHIKLMNPQRSTVWY][0-9]{1,4}[ACDEFGHIKLMNPQRSTVWY]$' # Search for single mutant variants, up to position length 9999
    prismdata.dataframe = prismdata.dataframe[prismdata.dataframe["variant"].map(lambda x: bool(re.match(r, x)))] # True if match, false else
    return(prismdata)


# Merge files
for mave_key in mave_path_dict:

    # Extract correct matching computational datasets for each MAVE experiment
    try:
        prismdata_other_paths = [path_dict_list[i][mave_key][0] for i in range(1, file_count)]
    except:
        prismdata_other_paths = ""
        print("Unable to find all non-MAVE data for:", mave_key)
        for line in [path_dict_list[i][mave_key] for i in range(1, file_count)]: print(line)
        continue

    for mave_path in mave_path_dict[mave_key]:
        outname = os.path.basename(mave_path)
        outname = outname.replace("prism_mave", "prism_merged")
        print("\n", outname)

        # Load prismdata mave
        prismdata_mave = prismparser.read(mave_path)
        print(prismdata_mave)

        # Remove multivariants
        prismdata_mave = remove_multivariants_prismdata(prismdata_mave)

        # Load Gemme, Rosetta, DSSP prismdata files, excluding prismdata mave
        prismdata_other_list = [prismparser.read(path) for path in prismdata_other_paths]

        # 4 quant normalize before merging
        from sklearn.preprocessing import QuantileTransformer
        if verbose: print("Rank normalizing GEMME only")

        for i in range(len(prismdata_other_list)):
            pdf = prismdata_other_list[i]
            if "gemme_score" in pdf.dataframe.columns:
                qt = QuantileTransformer(n_quantiles = len(pdf.dataframe))
                X_norm = qt.fit_transform(pdf.dataframe[["gemme_score"]])
                pdf.dataframe["gemme_score"] = X_norm
                if verbose:  print(X_norm)

        # Merge
        prismdata_merged = prismdata_mave.merge(prismdata_other_list, min_identity=0.0001, min_coverage=0.0001, mismatch="remove", allow_inserts=True, allow_deletions=True, merge = "outer")

        # Order
        prismdata_merged.order_variants(aa_order="ACDEFGHIKLMNPQRSTVWY=*~")

        # Save
        print("Saving to", outdir+"/"+outname, "\n")
        prismparser.write(outdir+"/"+outname, prismdata_merged)

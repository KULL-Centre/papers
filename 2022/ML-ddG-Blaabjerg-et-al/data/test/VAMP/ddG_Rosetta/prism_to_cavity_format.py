# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 10:20:21 2020

@author: Lasse
"""
import os
import pandas as pd
import PrismData

# Initialize
parser = PrismData.PrismParser()
df_out = pd.DataFrame(columns = ['pdbid', 'chainid', 'variant', 'score']) 

# Load data
file_list = [f for f in os.listdir(os.getcwd()) if f.endswith('.txt')]

# Read files
for file in range(len(file_list)):

    file_path = os.path.join(os.getcwd(), file_list[file])
    df_in = parser.read(file_path)
    df_tmp = pd.DataFrame(columns = ['pdbid', 'chainid', 'variant', 'score'])
    
    df_tmp['variant'] = df_in.dataframe["variant"]
    df_tmp['score'] = df_in.dataframe["Rosetta_ddg_score"]
    df_tmp['pdbid'] = df_in.metadata["protein"]["pdb"]
    df_tmp['chainid'] = df_in.metadata["protein"]["chain_id"]

    df_out = df_out.append(df_tmp)    

# Drop "=" variants
df_out = df_out[df_out["variant"].str[-1] != "="]

# Save
df_out.to_csv("ddg.csv", index=False)

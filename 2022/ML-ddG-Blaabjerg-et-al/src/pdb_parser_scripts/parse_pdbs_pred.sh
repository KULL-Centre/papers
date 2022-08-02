#!/bin/bash

# Settings
counter=1
dir=$(pwd)
reduce_exe=$dir/pdb_parser_scripts/reduce/reduce_src/reduce
pdb_dir=$1
pdbs=$pdb_dir/raw/*.pdb
n_pdbs=$(echo $pdbs | wc -w)

# Create data directories
mkdir -p $pdb_dir/cleaned
mkdir -p $pdb_dir/parsed

# Clean pdbs
for pdb in $pdbs;
do
    python $dir/pdb_parser_scripts/clean_pdb.py --pdb_file_in $pdb  \
                                                --out_dir $pdb_dir/cleaned/ \
                                                --reduce_exe $reduce_exe #&> /dev/null

    # Check for exit code 0 and skip file if not 0.
    if [ $? -eq 0 ]
    then
    echo "Successfully cleaned $pdb. $counter/$n_pdbs."
    else
    echo "Error when cleaning $pdb. Skipping.." >&2
    fi
    counter=$((counter+1))
done

# Parse pdbs and save in npz format
counter=1
out_dir=$pdb_dir/parsed
for pdb_clean in $pdb_dir/cleaned/*.pdb
do
    python $dir/pdb_parser_scripts/extract_environments.py --pdb_in $pdb_clean  \
                                                           --out_dir $out_dir  &> /dev/null

    # Check for exit code 0 and skip file if not 0.
    if [ $? -eq 0 ]
    then
    base_pdb=$(basename $pdb_clean)
    echo "Successfully parsed $pdb_clean. \
Finished $counter/$n_pdbs."
    else
    echo "Error extracting $pdb_clean. Skipping.." >&2
    fi
    counter=$((counter+1))
done

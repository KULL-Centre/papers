#!/bin/bash

# Settings
reduce_exe=reduce/reduce

# Clean pdbs
counter=1
dir=$(pwd)
parentdir=$(dirname $(pwd))
reduce_exe=$dir/pdb_parser_scripts/reduce/reduce
pdb_dir=$1
pdb=$2
pdbs=$parentdir/$pdb_dir/raw/*.pdb
n_pdbs=$(echo $pdbs | wc -w)

do
    python $dir/pdb_parser_scripts/clean_pdb.py --pdb_file_in $pdb  \
                                                --out_dir $parentdir/$pdb_dir/cleaned/ \
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
out_dir=$parentdir/$pdb_dir/parsed
pdb_clean="$pdb_clean"
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

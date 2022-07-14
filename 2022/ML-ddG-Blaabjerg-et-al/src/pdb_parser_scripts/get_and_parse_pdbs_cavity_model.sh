#!/bin/bash
conda_activate_path = $1
source $conda_activate_path cavity-model

# Settings
pdbid_filename = $2
reduce_exe = reduce/reduce
n_pdbs = $(cat $pdbid_filename | wc -l)

# Create data directories
mkdir -p data/cavity/structure/raw
mkdir -p data/cavity/structure/cleaned
mkdir -p data/cavity/structure/parsed

# Clean pdbid text file
tr -d '\n' < $pdbid_filename > pdbid_list_clean

# Download the protein structures (.pdb files) from the protein data bank
counter=1
while read p;
do
    # Download
    wget -O data/cavity/structure/raw/$p.pdb \
    'http://www.pdb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId='$p\
    &> /dev/null

    # Check that wget gives exit code 0. If not, remove potentially non-existant or partial PDB file.
    if [ $? -eq 0 ]
    then
        echo "Successfully downloaded $p.pdb."
    else
        echo "Could not download $p.pdb. Attempting to remove potentially non-existant or \
potentially partial $p.pdb" >&2
    
        rm data/cavity/structure/raw/$p.pdb 
    fi
    counter=$((counter+1))
done < pdbid_list_clean


# Clean all downloaded pdbs
counter=1
for pdb in data/cavity/structure/raw/*.pdb;
do
    python pdb_parser_scripts/clean_pdb.py --pdb_file_in $pdb  \
                                           --out_dir data/cavity/structure/cleaned/ \
                                           --reduce_exe $reduce_exe &> /dev/null

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
for pdb_clean in data/cavity/structure/cleaned/*.pdb
do
    python pdb_parser_scripts/extract_environments.py --pdb_in $pdb_clean --out_dir data/cavity/structure/parsed &> /dev/null

    # Check for exit code 0 and skip file if not 0.
    if [ $? -eq 0 ]
    then
    base_pdb=$(basename $pdb_clean)
    echo "Successfully parsed $base_pdb. \
Finished $counter/$n_pdbs."
    else
    echo "Error extracting $pdb_clean. Skipping.." >&2
    fi
    counter=$((counter+1))
done

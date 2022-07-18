import os
from Bio import PDB

pdb_io = PDB.PDBIO()
pdb_parser = PDB.PDBParser()

pdbfile_list = [f for f in os.listdir(os.getcwd()) if f.endswith(".pdb")]

for pdbfile in pdbfile_list:

    print(pdbfile)

    structure = pdb_parser.get_structure(" ", pdbfile)

    for model in structure: # Hack to change res numbers to very large values so there are no duplicate numbering later
        for chain in model:
            for i, residue in enumerate(chain.get_residues()):
                res_id = list(residue.id)
                res_id[1] = int(i+10000)
                residue.id = tuple(res_id)

    for model in structure:
        for chain in model:
            for i, residue in enumerate(chain.get_residues()):
                res_id = list(residue.id)
                res_id[1] = int(i+1)
                residue.id = tuple(res_id)


    pdb_io.set_structure(structure)
    pdb_io.save(pdbfile)

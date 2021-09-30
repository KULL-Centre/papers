Here is an example of how to rescale protein-water interactions in Martini 3 using the provided python script.

You will need a topology file (.top) for your system. For this example the provided topology file is A1_system_topology.top.

Your protein must be named something with "Protein" under [ moleculetype ] in the topology, e.g. the protein A1 in this example is named "Protein 1" in this example.

You can run the script from the command line using:
```
python PW_rescaling_martini3.py -i A1_system_topology.top -o A1_system_topology_lambda1.06.top -l 1.06
```
This generates a new topology file A1_system_topology_lambda1.06.top as output, where Lennard-Jones interactions between protein and water beads have been rescaled by a factor 1.06.

WARNING! Be careful if you have other molecules than proteins (e.g. lipids) in your system, as any bead-types that are also present in the protein will have rescaled interactions with water. You could potentially get around this by renaming the bead-types of the protein.

Also, if you have multiple different proteins, you must specify the number of proteins with the flag -n. Make sure they are all named something with "Protein" under [ moleculetype ].

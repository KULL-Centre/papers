#PyMOL surface area batch script
#Frederik Theisen, 2020
#Usage: place file in active folder and execute "run <filename.extension>" #Example: run script.txt
#Ensure the existence of properly formatted database file “[path/]database.txt”
#Database format: <pdb id> <idp chain id>
#Database format example: “1AXC	B”
#Set active folder: cd <folder path>

print("PyMOL IDP COMPLEX SURFACE AREA SCRIPT")

DATABASE = 'id_complex_database.txt'	#Database path (path relative to working dir)
MINIMUMLENGTH = 4 			#Minimum IDP chain length
QUALITY = 2				#Surface calculation resolution (0 = very fast, 4 = very slow, higher accuracy. Default = 2) 

#Define function to get surface area of structure
def SASA():
	#Set variables (1.4 Å is standard water molecule size)
	cmd.do("set solvent_radius, 1.4")
	cmd.do("set dot_density, " + str(QUALITY))

	#Prepare molecule
	cmd.do("remove hydro") 		#remove hydrogens
	cmd.do("remove solvent") 	#remove water solvent
	cmd.do("select s, hetatm") 	#select non-protein atoms
	cmd.do("remove s") 		#remove non-protein selection
	cmd.do("show dots") 		#prepare protein calculation surface 
	cmd.do("set dot_solvent, on")	#approximate solvent size

	
	#Count models
	states = cmd.count_states("all")

	total = 0
	polar = 0
	nonpolar = 0

	#Loop through states and calculate surface area
	for i in range(0,states):
		total += cmd.get_area("all", i+1)
		polar += cmd.get_area("elem N+O", i+1)
		nonpolar += cmd.get_area("elem C+S", i+1)

	#Find average surface areas
	avg_total = total/states
	avg_polar = polar/states
	avg_nopol = nonpolar/states

	return [avg_total, avg_polar, avg_nopol]

print("Starting Script")
print("Loading database...")

with open(DATABASE) as f:
	lines = f.readlines()
	
print("Database loaded")
print("Starting calculations...")
print("")

#Output container
output = []

for i in range(0, len(lines)):
	#Clear PyMOL
	cmd.delete("all")
	
	#Get line
	line = lines[i].strip()

	#Get PDB code
	pdbid = line[:4].strip()

	#Fetch molecule
	cmd.fetch(pdbid)

	#Get IDP chain identifier
	disorderedchainref = line[5:].strip()
	
	#Select IDP with only visible amino acids and save as "idp"
	cmd.do("select idp, chain " + disorderedchainref + ", 0, 1, 0, 1")
	#Get IDP sequence and remove chain ID string (selection string, state = only visible, quiet = yes)
	sequence = cmd.get_fastastr("idp")[5:].strip().replace("\n","")

	#Skip the complex if the IDP contains modified or non-standard residues
	if '?' in sequence: continue

	#Filters short IDP chains
	if len(sequence) < MINIMUMLENGTH: continue

	#Get complex surface (SASA defined in function, returns array containing total, polar and nonpolar surface area)
	surface_complex = SASA()

	#Remove IDP
	cmd.do("remove chain " + str(disorderedchainref))

	#Get folded protein surface
	surface_folded = SASA()

	#Converts values to strings and saves in output container
	output.append([pdbid, sequence, str(surface_folded[0]), str(surface_folded[1]), str(surface_folded[2]), str(surface_complex[0]), str(surface_complex[1]), str(surface_complex[2])])

print("")
print("Calculations completed")
print("Results: " + str(len(output)))
print("Results format: PDB, IDP seq, folded domain (fd) total surface, fd polar, fd non-polar, complex total, complex polar, complex non-polar")
print("Saving results to file...")

f = open('output.txt', 'w+')

for e in output:
	f.write(str(e[0]) + "\t" 
	+ str(e[1]) + "\t" 
	+ str(e[2]) + "\t" 
	+ str(e[3]) + "\t" 
	+ str(e[4]) + "\t" 
	+ str(e[5]) + "\t" 
	+ str(e[6]) + "\t" 
	+ str(e[7]) + "\n")

f.close()

print("Results saved in file: output.txt")

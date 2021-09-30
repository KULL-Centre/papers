import sys
import numpy as np

#Topology file with rubberbands
topfile="all_PRO.top"

#Output topology with modified rubberbands
outfile="rubberbands_all_PRO.top"

#Number of domains in which to keep interdomain rubberbands
nr_domains = 2

#List of start and stop atoms of domain borders
#Given in the form [domain 1 start-atom, domain 1 stop-atom, domain 2 start-atom, domain 2 stop-atom]
start_stop_atoms = [23, 207, 244, 425]

#Check that there are start and stop atoms for all domains
assert len(start_stop_atoms) == 2*nr_domains

#Function to return rubberbands only if both atoms are within atom_start and atom_stop
#Toplines are all lines of the topology file - start_line and end_line are the indeces of the rubber-band section
def rubberbands_to_keep_domain(atom_start, atom_stop, toplines, start_line, end_line):
    
    lines_to_keep = []

    for line in toplines[start_line:end_line]:
        linesplit = line.split()

        if int(linesplit[0]) >= int(atom_start) and int(linesplit[0]) <= int(atom_stop) and int(linesplit[1]) >= int(atom_start) and int(linesplit[1]) <= int(atom_stop):
            lines_to_keep.append(line)
    return lines_to_keep

#Find Rubber band lines in topfile
with open(topfile, 'r') as f:
    toplines = f.readlines()

for i in range(len(toplines)):
    if '; Rubber band' in toplines[i]:
        start_line = int(i+1)
        break

for i in range(len(toplines)):
    if ';' in toplines[i] and i>start_line:
        end_line = int(i-1)
        break

print("Rubber bands are from line %s to line %s in %s" % (start_line, end_line, topfile))


#Write to outfile
with open(outfile, 'w') as f:
    
    #Write topology before Rubberbands
    for line in toplines[:start_line]:
        f.write(line)
    
    #Get rubberbands for each domain and write to file
    j=0
    for i in range(nr_domains):
        print("Writing Rubber bands for domain %i from atom %i to %i" % (i+1, start_stop_atoms[j], start_stop_atoms[j+1]))
        rubberbands = rubberbands_to_keep_domain(start_stop_atoms[j], start_stop_atoms[j+1], toplines, start_line, end_line)
        for line in rubberbands:
            f.write(line)
        
        j+=2
    
    #Write rest of topology
    for line in toplines[end_line:]:
        f.write(line)

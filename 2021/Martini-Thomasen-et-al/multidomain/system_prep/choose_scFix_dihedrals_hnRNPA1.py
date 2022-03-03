topfile="rubberbands_all_PRO.top"
outfile="dihedrals_rubberbands_all_PRO.top"

nr_domains = 2
start_stop_atoms = [23, 207, 244, 425]
assert len(start_stop_atoms) == 2*nr_domains

with open(topfile, 'r') as f:
    toplines = f.readlines()

def dihedrals_to_keep_domain(atom_start, atom_stop, toplines, start_line, end_line):

    lines_to_keep = []

    for line in toplines[start_line:end_line]:
        linesplit = line.split()

        if int(linesplit[0]) >= int(atom_start) and int(linesplit[0]) <= int(atom_stop) and int(linesplit[1]) >= int(atom_start) and int(linesplit[1]) <= int(atom_stop) and int(linesplit[2]) >= int(atom_start) and int(linesplit[2]) <= int(atom_stop):
            lines_to_keep.append(line)
    return lines_to_keep

for i in range(len(toplines)):
    if '; SC-BB-BB and BB-BB-SC scFix' in toplines[i]:
        start_line_SCBBBB_BBBBSC = int(i+1)
        break
        
for i in range(len(toplines)):
    if ';' in toplines[i] and 'BB-BB' not in toplines[i] and i>start_line_SCBBBB_BBBBSC:
        end_line_SCBBBB_BBBBSC = int(i-1)
        break
        
print("SC-BB-BB and BB-BB-SC scFix are from line %s to line %s in %s" % (start_line_SCBBBB_BBBBSC, end_line_SCBBBB_BBBBSC, topfile))

for i in range(len(toplines)):
    if '; SC-BB-BB-SC scFix' in toplines[i]:
        start_line_SCBBBBSC = int(i+1)
        break
        
for i in range(len(toplines)):
    if '[' in toplines[i] and i>start_line_SCBBBBSC and 'BB-BB' not in toplines[i]:
        end_line_SCBBBBSC = int(i-1)
        break
        
print("SC-BB-BB-SC scFix are from line %s to line %s in %s" % (start_line_SCBBBBSC, end_line_SCBBBBSC, topfile))

with open(outfile, 'w') as f:

    #Write topology before scFix dihedrals
    for line in toplines[:start_line_SCBBBB_BBBBSC]:
        f.write(line)

    j=0
    for i in range(nr_domains):
        dihedrals = dihedrals_to_keep_domain(start_stop_atoms[j], start_stop_atoms[j+1], toplines, start_line_SCBBBB_BBBBSC, end_line_SCBBBB_BBBBSC)
        for line in dihedrals:
            f.write(line)

        j+=2

    #Write topology from first scFix dihedrals to next scFix dihedrals 
    for line in toplines[end_line_SCBBBB_BBBBSC:start_line_SCBBBBSC]:
        f.write(line)
    
    j=0
    for i in range(nr_domains):
        dihedrals = dihedrals_to_keep_domain(start_stop_atoms[j], start_stop_atoms[j+1], toplines, start_line_SCBBBBSC, end_line_SCBBBBSC)
        for line in dihedrals:
            f.write(line)
        j+=2

    #Write rest of topology
    for line in toplines[end_line_SCBBBBSC:]:
        f.write(line)



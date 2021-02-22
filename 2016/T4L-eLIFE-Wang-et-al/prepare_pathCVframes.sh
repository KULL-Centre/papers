#!/bin/sh
#

for i in `seq 1 1 31`
do
        ipdb=`printf %03d $i`

cat <<EOF >pdb2frame.tcl
mol new pdb/${ipdb}.pdb type pdb first 0 last -1 step 1 waitfor 1
set all [atomselect top all]
\$all set occupancy 0
\$all set beta 0
#set align [atomselect top "backbone"]
set align [atomselect top "name CA"]
\$align set occupancy 1.0
set bias [atomselect top "resid 99 to 126 and not hydrogen"]
\$bias set beta 1.0
#set all [atomselect top "name CA or (resid 99 to 126 and not hydrogen)"]
set all [atomselect top "(name CA and resid 95 to 126) or (resid 99 to 126 and not hydrogen)"]
\$all writepdb tmp.pdb
quit
EOF
vmd -dispdev text -e pdb2frame.tcl 2>&1 >/dev/null
echo /usr/local/vmd-1.9.1/bin/vmd -dispdev text -e pdb2frame.tcl
perl Moil2Charmm.pl -G tmp.pdb -F T4L_charmm.pdb -O frame${i}.pdb
done

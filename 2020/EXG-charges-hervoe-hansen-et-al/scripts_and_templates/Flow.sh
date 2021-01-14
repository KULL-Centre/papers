#!/bin/sh
# Stefan H. Hansen & Yong Wang

##############
# FLOW CHART #
##############

# Source your amber.sh where AMBER is installed.
source /temp1/Software/amber16/amber.sh

# Step 1
# Modify the residue names of the residues you wish to titrate
# Note that Amber currently does not allow to titrate ARG, N-, and C- terminal residues
perl -pi -e's/HIS/HIP/,s/ASP/AS4/,s/GLU/GL4/,s/CYS/CYX/' EXG_WT_Fixed.pdb 

# Step 2
# Prepare the topology file (prmtop) and input coordinate (inpcrd) file using tleap. 
cat <<EOF >tleap.in
source leaprc.constph
EXGwt = loadPDB ../PDB_files/EXG_WT_Fixed.pdb
loadAmberParams frcmod.ionsjc_tip3p

# Make the disulfide bonds
bond EXGwt.8.SG EXGwt.107.SG
solvateOct  EXGwt TIP3PBOX 10.0
addIonsRand EXGwt Cl- 3 Na+ 1

saveAmberParm EXGwt EXGwt.parm7 EXGwt.rst7
quit
EOF
tleap -f tleap.in

IonCon=`echo 3 | awk '{printf("%12.8f \n",$1/(0.6*6.4^3))}'`

# Step 3
# Creating the CPIN file
echo /usr/bin/python /sbinlab/hervoe/software/AMBER/amber14/bin/cpinutil.py --help
/usr/bin/python /sbinlab/hervoe/Software/AMBER/amber14/bin/cpinutil.py -resnames HIP AS4 GL4 LYS TYR -p EXGwt.parm7 -o EXGwt.cpin -op EXGwt.explicit.parm7
exit

# Step 4
# Minimization (min.mdin & min.sh [GPU])
mpirun -n 28 pmemd.MPI -O -i ../scripts_and_templeates/min.mdin -p EXG.explicit.parm7 -c EXG.rst7 -o EXG.min.mdout -r EXG.min.rst7 -ref EXG.rst7 -cpin EXG.cpin # CPU

# Step 5
# Heating (heat.mdin & heat.sh [GPU])
mpirun -n 28 pmemd.MPI -O -i ../scripts_and_templeates/heat.mdin -c EXG.min.rst7 -p EXG.explicit.parm7 -ref EXG.min.rst7 -cpin EXG.cpin -o EXG.heat.mdout -r EXG.heat.rst7 -x EXG.heat.nc #CPU

# Step 6
# Equilibration (Relaxation) (equil.mdin & equil.sh [GPU])
mpirun -n 28 pmemd.MPI -O -i ../scripts_and_templeates/equil.mdin -p EXG.explicit.parm7 -c EXG.heat.rst7 -ref EXG.min.rst7 -cpin EXG.cpin -o EXG.equil.mdout -cpout EXG.equil.cpout -r EXG.equil.rst7 -x EXG.equil.nc -cprestrt EXG.equil.cpin #CPU

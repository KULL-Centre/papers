import os
import mdtraj as md

cg_traj = '../prodrun_nopbc.xtc'
cg_top = '../../../two_FUS_init/PRO_CG.gro'

BM='/home/projects/ku_10001/people/fretho/software/backward-v5/initram-v5_200.sh'
gmx='/home/projects/ku_10001/apps/gromacs-2018.1_AVX2_256_12JUN2018/bin/gmx_mpi'

outfile_xtc = '../prodrun_AAbackmapped.xtc'
outfile_gro = '../prodrun_AAbackmapped.gro'

#Load CG trajectory
traj = md.load(cg_traj, top=cg_top)

#Loop over CG trajectory frames
#for i in range(0,len(traj),4):

for i in range(0,len(traj)):
    
    while os.path.isfile(f'AA_frame{i}.pdb') == False:
         
        #Save frame to gro file
        md.Trajectory.save_gro(traj[i], 'CG_frame.gro')

        #Backmap
        os.system('bash ' + BM + ' -f CG_frame.gro -p topol.top -to charmm36 >/dev/null')

        os.system(str(gmx) + ' editconf -f backmapped.gro -o AA_frame%s.pdb -quiet' % str(i))
        #Remove files for next iteration
        os.system('rm backmapped.gro')
        os.system('rm CG_frame.gro')
	
#Load initial AA frame as traj
AA_traj = md.load('AA_frame0.pdb')

#Loop over all frames but first and concatenate frame to traj
for i in range(0,len(traj)):
    frame = md.load('AA_frame%i.pdb' % i)
    frame.time = i*1000
    AA_traj = md.Trajectory.join(AA_traj, frame)

#Center protein
AA_traj_centered = md.Trajectory.center_coordinates(AA_traj)

#Save xtc and gro files
md.Trajectory.save_xtc(AA_traj_centered, outfile_xtc)
md.Trajectory.save_gro(AA_traj_centered[0], outfile_gro)

print("Some info on the generated trajectory: " + str(AA_traj_centered))

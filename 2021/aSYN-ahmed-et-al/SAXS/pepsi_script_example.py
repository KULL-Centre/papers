import subprocess
import os
import sys
import numpy as np

import os

Pepsi_path = '/lustre/hpc/sbinlab/software/Pepsi-SAXS-SANS/Pepsi-SAXS'
SAXS_data = 'exp.dat'

nr_of_frames = 0
for _, _, files in os.walk('.'):
        for file in files:
                if file.startswith('frame') and file.endswith('.pdb'):
                        nr_of_frames += 1
                        
#first, fit all frames and all parameters
step=1
if step==1:
	for i in range(0,nr_of_frames+1):
		filename = "frame_%d.pdb" % i
		outfile = "frame_%d.fit" % i
		subprocess.run([Pepsi_path, filename, SAXS_data , '-o '+  outfile, '-ms ' +' 2' ,'-cst' ,'-json'])

#find default r0
subprocess.run([Pepsi_path, 'frame_10.pdb', SAXS_data , '-o find_r0.fit', '-ms ' +' 2', '--r0_min_factor 1.0', '--r0_max_factor 1.0', '--r0_N 1', '-cst', '-json'])
with open('find_r0.log', 'r') as f:
        lines = f.readlines()
co=0
for line in lines:
        if "Best r0 found" in line:
                co+=1 #1
        if co==1 and "value" in line:
                co+=1 #2
                default_r0=float(line.split('"')[3])


#then, parse all the log files to find min/max r0 factors and best d_rho
print('parsing', nr_of_frames+1, 'log files')
r0_list = []
d_rho_list = []
#print (d_rho_list)
c=0
for i in range(0,nr_of_frames+1):
        filename = "frame_%d.log" % i
        try:
                with open(filename, 'r') as f:
                        lines = f.readlines()
                for line in lines:
                        if "Best d_rho found" in line:
                                c+=1 #1
                        if c==1 and "value" in line:
                                c+=1 #2
                                value=float(line.split('"')[3])
                                d_rho_list.append(value)
                        if "Best r0 found" in line:
                                c+=1 #3
                        if c==3 and "value" in line:
                                c+=1 #4
                                value=float(line.split('"')[3])/default_r0
                                r0_list.append(value)
        except:
                print(filename, 'not found, skipping it.')
                
                               
                

#Take the average to find ensemble min/max r0 factors and d_rho
r0 = str(np.average(r0_list))
d_rho = str(np.average(d_rho_list)*10)



#refit all of the frames, using the ensemble parameters
for i in range(0,nr_of_frames+1):
        filename = "frame_%d.pdb" % i
        outfile = "Ave_AA_frame%d.fit" % i
        subprocess.run([Pepsi_path, filename, SAXS_data , '-o ' + outfile, '--dro ' + d_rho, '-ms ' +' 2', '--r0_min_factor '+ r0, '--r0_max_factor '+ r0, '--r0_N 1', '-cst'])
        
print('ensemble d_rho', d_rho)
print('ensemble r0', r0)
print('default r0', default_r0)
print('done')

# import the reweighting script
import matplotlib.pyplot as plt
import numpy as np
import bme_reweight as bme
import sys

# chooce contrast with j. 0:SAXS, 1:SANS 0%, 2:SANS 42%, 3: SANS 70%
try:
    j = int(sys.argv[1])
except:
    j = 0

#define theta list
if j == 0:
    t_list = [50,100,200,300,400,500,1000,2000,3000]
elif j == 1:
    t_list = [1,10,50]
elif j == 2:
    t_list = [1,10,50]
elif j == 3:
    t_list = [1,10,50]
elif j > 3:
    t_list = [1,10,50]

# define name of experimental datafiles
exp_file = 'exp_data_%d.dat' % j

# define name of experimental quantities # calculated from the trajectory
sim_file = 'calc_data_%d.dat' % j

# initialize reweighting class 
rew = bme.Reweight()

# load data
rew.load(exp_file,sim_file)

# do optimization
Neff_list = []
chia_list= []
chib_list = []
weights = {}
for t in t_list:
    outfile1 = 'BME_sas_theta%d_%d' % (t,j)       
    chib,chia, srel = rew.optimize(theta=t) # chib = chi before optimization, chia = chi after optimization
    chi2_b,chi2_a = rew.weight_exp(exp_file,sim_file, outfile1) #save fitting data for each data
    w_opt = rew.get_weights()
    weights[t]=w_opt
    chia_list.append(chia)
    chib_list.append(chib)
    Neff_list.append(np.exp(srel))
data = t_list, Neff_list, chib_list, chia_list

#save weights for each theta
best_theta = t_list[0]
for t in t_list:
    outfile2 = 'weights_theta%d_%d.dat' % (t,j)
    to_write=[]
    for i in range(0, len(weights[t])):
            to_write.append('frame_%d\t\t%f\n' % (i, weights[t][i]))
    with open(outfile2, 'w') as f:
       for w in to_write:
            f.write(w)
    
# save data
Ni = len(Neff_list)
filename = 'BME_%d.dat' % j
with open(filename, 'w') as f:
    f.write("Theta    \t Neff    \t Chi2r_before \t Chi2r_after \n")
    for i in range(Ni):
        f.write("%f \t %f \t %f \t %f \n" % (t_list[i],Neff_list[i],chib_list[i],chia_list[i]))

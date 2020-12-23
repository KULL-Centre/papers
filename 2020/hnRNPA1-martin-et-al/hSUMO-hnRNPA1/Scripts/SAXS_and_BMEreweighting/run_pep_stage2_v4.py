import subprocess
import os
import sys
import numpy as np
from scipy.optimize import curve_fit

# user input
first_frame = 1000
last_frame = 20000
first_contrast = 0
last_contrast = 0
first_combination = 1 # for BME
last_combination = 4 # for BME
skip_last = 0 # skip last points (from SAXS dataset)

# define paths
PepsiSAXS_path = '/lindorffgrp-isilon/andreas/programs/Pepsi-SAXS_linux_2.4 '
PepsiSANS_path = '/lindorffgrp-isilon/andreas/programs/Pepsi-SANS_linux_2.4 '
SAXS_data = '/lindorffgrp-isilon/thomasen/hnRNPA1/dHLCD_SAXS_data/FLdH_50mMHEPES_%imMNaCl_reduced' % str(sys.argv[1])
#SANS0_data = 'TIA-21-0'
#SANS42_data = 'TIA-21-42'
#SANS70_data = 'TIA-21-70'

# print
nr_of_frames = last_frame-first_frame+1
print('\n\n nr_of_frames = %d' % nr_of_frames)

for j in range(first_contrast,last_contrast+1):
    
    # set flags depending on contrast
    flags = ' -cst'
    if j == 0:
        Pepsi_path = PepsiSAXS_path
        Data = SAXS_data
        extra_name = ''
        Ifit_column = 3
    else:
        Pepsi_path = PepsiSANS_path
        flags += ' --deuterated A --deut 1.00'
        Ifit_column = 4
        if j == 1:
            Data = SANS0_data
            extra_name = '_deut100_d2o0'
            flags += ' --d2o 0.0'
        if j == 2:
            Data = SANS42_data
            extra_name = '_deut100_d2o42'
            flags += ' --d2o 0.42'
        if j == 3:
            Data = SANS70_data
            extra_name = '_deut100_d2o70'
            flags += ' --d2o 0.7'
    Data_path = '/lindorffgrp-isilon/thomasen/hnRNPA1/dHLCD_SAXS_data/' + Data + '.dat'
    
    # import parameters from stage 1
    params_file = 'params_' + str(j)
    r0_0,drho_0,I0_0,cst_0 = np.genfromtxt(params_file,skip_header=1,skip_footer=0,usecols=[0,1,2,3], unpack=True)

    # import thetas
    bme_file = 'BME_%d.dat' % j
    thetas = np.genfromtxt(bme_file,skip_header=1,skip_footer=0,usecols=[0],unpack=True)

    # initiate summary file
    summary_file = 'summary_file_%d.dat' % j
    with open(summary_file,'w') as f:
        f.write('theta       chi2r         chi2r_w           chi2r_w_opt         r0          drho          I0           cst          S_opt         B_opt\n')
    
    # calculate scattering curves for each theta value (from BME) and find minimum of chi2r_w
    chi2r_w_opt_min = 999999 # a very large number
    for theta in thetas:
        # import weights from BME
        weights_file = 'weights_theta%d_%d.dat' % (theta,j)
        w = np.genfromtxt(weights_file,skip_header=0,skip_footer=0,usecols=[1], unpack=True)
        
        # calculate weighted averages of parameters
        r0=str(np.sum(w*r0_0))
        drho=str(np.sum(w*drho_0))
        I0=str(np.sum(w*I0_0))
        cst=str(np.sum(w*cst_0))
        
        Ifit_sum = 0.0
        Ifit_w   = 0.0

        # fit all frames with fixed parameters
        for i in range(first_frame,last_frame+1):
            filename = "../Backmapping_frames/AA_frame%d.pdb" % i
            outfile = "w_AA_frame%dchain_%d.fit" % (i,j)
            extra_flags = ' -o ' + outfile + ' --dro ' + drho + ' --r0_min_factor ' + r0 + ' --r0_max_factor ' + r0 + ' --r0_N 1' + ' --cstFactor ' + cst + ' --I0 ' + I0
            command = Pepsi_path + filename + " " + Data_path + flags + extra_flags
            os.system(command)
            q,Iexp,dI,Ifit =  np.genfromtxt(outfile,skip_header=6,skip_footer=0,usecols=[0,1,2,Ifit_column],unpack=True)
            Ifit_sum += Ifit
            Ifit_w   += w[i]*Ifit

        # find optimal I0 and cst
        def func(q,S,B):
            I = S*Ifit_w+B
            return I
        S0 = 1.00
        B0 = 0.00
        popt,pcov = curve_fit(func,q,Iexp,sigma=dI,p0=[S0,B0])
        S_opt = popt[0]
        B_opt = popt[1]

        # calculate chi2r and weighted chi2r
        Ifit_avg = Ifit_sum/nr_of_frames
        Ifit_w_opt = S_opt * Ifit_w + B_opt
        R = (Ifit_avg - Iexp)/dI
        R_w = (Ifit_w - Iexp)/dI
        R_w_opt = (Ifit_w_opt - Iexp)/dI
        K = 4 # no of parameters fitted
        M = len(Iexp) # nu of data points
        nu = M-K # degrees of freedom
        chi2r = sum(R**2)/nu
        chi2r_w = sum(R_w**2)/nu
        chi2r_w_opt = sum(R_w_opt**2)/nu

        # write to summary file
        with open(summary_file,'a') as f:
            f.write('%f     %f     %f      %f      %s     %s     %s     %s    %1.4f    %1.4e\n' % (theta,chi2r,chi2r_w,chi2r_w_opt,r0,drho,I0,cst,S_opt,B_opt))

        # find min chi2r_w and corresponding theta
        if chi2r_w_opt < chi2r_w_opt_min:
            chi2r_w_opt_min = chi2r_w_opt
            chi2r_w_min = chi2r_w
            chi2r_min = chi2r
            r0_min = r0
            drho_min = drho
            I0_min = I0
            cst_min = cst
            theta_min = theta
            S_opt_min = S_opt
            B_opt_min = B_opt
         
    # write to summary file
    with open(summary_file,'a') as f:
        f.write('%f     %f     %f     %f      %s     %s     %s     %s   %1.4f    %1.4e\n' % (theta_min,chi2r_min,chi2r_w_min,chi2r_w_opt_min,r0_min,drho_min,I0_min,cst_min,S_opt_min,B_opt_min))

    # repeat calculation for theta at minimum chi2r_w, and with refitted I0 and cst
    weights_file = 'weights_theta%d_%d.dat' % (theta_min,j)
    w = np.genfromtxt(weights_file,skip_header=0,skip_footer=0,usecols=[1],unpack=True)
    r0=str(np.sum(w*r0_0))
    drho=str(np.sum(w*drho_0))
    I0=str(np.sum(w*I0_0)*S_opt_min)
    cst=str(np.sum(w*cst_0)-B_opt_min)
    for i in range(first_frame,last_frame+1):
        filename = "../Backmapping_frames/AA_frame%d.pdb" % i
        outfile = "w_AA_frame%dchain_%d.fit" % (i,j)
        extra_flags = ' -o ' + outfile + ' --dro ' + drho + ' --r0_min_factor ' + r0 + ' --r0_max_factor ' + r0 + ' --r0_N 1' + ' --cstFactor ' + cst + ' --I0 ' + I0
        command = Pepsi_path + filename + " " + Data_path + flags + extra_flags
        os.system(command)
    
print('\n\n Done with Pepsi part! \n')

### generate calc_data and exp_data files for BME

for combination in range(first_combination,last_combination+1):
 
    if combination == 0:
        D = [-1]
    if combination == 1:
        D = [0]
    if combination == 2:
        D = [42]
    if combination == 3:
        D = [70]
    if combination == 4:
        D = [-1,0,42,70]
    if combination == 5:
        D = [-1,0]
    if combination == 6:
        D = [-1,42]
    if combination == 7:
        D = [-1,70]
    if combination == 8:
        D = [0,42]
    if combination == 9:
        D = [0,70]
    if combination == 10:
        D = [42,70]
    if combination == 11:
        D = [-1,0,42]     
    if combination == 12:
        D = [-1,0,70] 
    if combination == 13:
        D = [-1,42,70]
    if combination == 14:
        D = [0,42,70]

    # define names of files
    outfile_exp = 'trash_exp_data_' + str(combination) + '.dat'
    outfile_calc = 'calc_data_stage2_' + str(combination) + '.dat'

    # print header to exp file
    with open(outfile_exp, 'w') as f:
        f.write("# DATA=SAXS PRIOR=GAUSS\n")

    # make first column of header to calc file
    header = "# label"

    for j in range(len(D)):
  
        #make EXP file
        footerlines = 0
        if D[j] == -1:
            data_filename = "/lindorffgrp-isilon/thomasen/hnRNPA1/dHLCD_SAXS_data/FLdH_50mMHEPES_%imMNaCl_reduced.dat" % str(sys.argv[1]) 
            headerlines = 0
        else:
            data_filename = "" % D[j]
            headerlines = 1
        footerlines = 0
        q,Iexp,dI = np.genfromtxt(data_filename,skip_header=headerlines,skip_footer=footerlines,usecols=[0,1,2], unpack=True)
        if D[j] == -1:
            q = q[0:-skip_last]
            Iexp = Iexp[0:-skip_last]
            dI = dI[0:-skip_last]

        with open(outfile_exp, 'a') as f:
            for i in range(0,len(q)):
                f.write("%e \t %e \t %e \n" % (q[i],Iexp[i],dI[i]))

        # write CALC_FILE header
        for q_value in q:
            header += " \t %e" % q_value

    # print header to calc file
    header += " \n"
    with open(outfile_calc,'w') as f:
        f.write(header)

    # print rest of calc file
    for i in range(first_frame,last_frame+1):
        frame_number = i+1 # index with 1 instead of 0 for frame numbers
        frame_line = "frame_%d" % frame_number
        for j in range(len(D)):
            headerlines = 6
            footerlines = 0
            Ifit_col = 4
            if D[j] == -1:
                contrast = 0
                Ifit_col = 3
            elif D[j] == 0:
                contrast = 1
            elif D[j] == 42:
                contrast = 2
            elif D[j] == 70:
                contrast = 3
            fit_filename = "w_AA_frame%dchain_%d.fit" % (i,contrast)
            q,Iexp,dIexp,Ifit = np.genfromtxt(fit_filename,skip_header=headerlines,skip_footer=footerlines,usecols=[0,1,2,Ifit_col], unpack=True)
            for Ifit_value in Ifit:
                frame_line += " \t %e" % Ifit_value
        frame_line += '\n'  
        with open(outfile_calc,'a') as f:             
            f.write(frame_line)
    
print('\n Done with make_calc_exp part! \n')

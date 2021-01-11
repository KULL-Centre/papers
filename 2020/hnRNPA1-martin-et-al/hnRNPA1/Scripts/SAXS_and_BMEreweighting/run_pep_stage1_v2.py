'''
This program uses pepsi SAXS/SANS to:
1) calculate scattering with free parameter values (for each frame)
2) calculate scattering with mean parameter values (for all frames)
3) calculate scattering with mean values for r0, drho, cst, but optimise for I0
4) make files for BME
'''

import subprocess
import os
import sys
import numpy as np
from scipy.optimize import curve_fit

# user input
first_frame = 0 # exclude first 0 frames (too affected by initial structure)
last_frame = 20000
first_contrast = 0
last_contrast = 0
first_combination = 0 # for BME
last_combination = 0 # for BME

# define paths
PepsiSAXS_path = '/lindorffgrp-isilon/andreas/programs/Pepsi-SAXS_linux_2.4 '
PepsiSANS_path = '/lindorffgrp-isilon/andreas/programs/Pepsi-SANS_linux_2.4 '
SAXS_data = '/lindorffgrp-isilon/thomasen/hnRNPA1/FLA1_SAXSdata_08042020_averaged_newest/' + str(sys.argv[1])
#SANS0_data = 'TIA-21-0'
#SANS42_data = 'TIA-21-42'
#SANS70_data = 'TIA-21-70'

# print number of frames
nr_of_frames = last_frame - first_frame + 1
print('\n\n nr_of_frames = %d' % nr_of_frames)

# define function for getting chi2r
def get_chi2r(FREE,Pepsi_path,Data_path,flags,Ifit_column,first_frame,last_frame,r0,drho,I0,cst):
    '''
    This sub function takes pepsi settings and params, and returns a chi2r

    ## input:
    FREE: fit all params freely
    Pepsi_path: path to pepsi program
    Data_path: path to data
    flags: flags used in pepsi
    Ifit_column: column for fit (3 for SAXS, 4 for SANS)
    first_frame: first frame to include
    last_frame: last frame to include 
    r0, drho, I0, cst: values of fitting parameters (strings)
    '''
    K = 4 # no of parameters fitted, for calculating chi2r
    Ifit_sum   = 0.0
    nr_of_frames = last_frame-first_frame+1
    for i in range(first_frame,last_frame+1):
        filename = "../Backmapping_frames/AA_frame%d.pdb" % i
        if FREE:
            outfile = "Free_AA_frame%dchain_%d.fit" % (i,j)
            extra_flags = ' -o ' + outfile
        else:
            outfile = "Fix_AA_frame%dchain_%d.fit" % (i,j)
            extra_flags = ' -o ' + outfile + ' --dro ' + drho + ' --r0_min_factor ' + r0 + ' --r0_max_factor ' + r0 + ' --r0_N 1' + ' --cstFactor ' + cst + ' --I0 ' + I0
        command = Pepsi_path + filename + " " + Data_path + flags + extra_flags
        os.system(command)
        q,Iexp,dI,Ifit =  np.genfromtxt(outfile,skip_header=6,skip_footer=0,usecols=[0,1,2,Ifit_column],unpack=True)
        Ifit_sum += Ifit
    Ifit_avg = Ifit_sum/nr_of_frames
    M = len(Iexp) # number of data points
    nu = M-K # degrees of freedom
    chi2r = sum(((Ifit_avg-Iexp)/dI)**2)/nu
    return chi2r,Ifit_avg,q,Iexp,dI

# write summary file
summary_file = 'summary_0_file.dat'
with open(summary_file,'w') as f:
    f.write('contrast     \tchi2r_mean     \tchi2r_opt      \tr0               \tdrho \tI0      \tcst     \tS_opt    \tB_opt\n')

# loop over contrasts
for j in range(first_contrast,last_contrast+1):

    # set settings and flags depending on contrast
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
    Data_path = Data + '.dat'
    
    # generate fit file
    fit_file = 'fit_file_%d.dat' % j
    with open(fit_file,'w') as f:
        f.write('#  Fitting with free parameter values:\n')
        f.write('#  I0_free        \tcst             \tchi2r \n')
   
    # fit all frames and all parameters
    chi2r_free = get_chi2r(1,Pepsi_path,Data_path,flags,Ifit_column,first_frame,last_frame,'dummy','dummy','dummy','dummy')[0]

    # write to fit file
    with open(fit_file,'a') as f:
        f.write('0.00000000        \t0.00000000000          \t%1.2f \n' % chi2r_free)
        
    # find default dro and r0
    extra_flags = ' -o find_r0_dro.fit --dro 1.0 --r0_min_factor 1.0 --r0_max_factor 1.0 --r0_N 1'
    filename = "../Backmapping_frames/AA_frame%d.pdb" % first_frame
    command = Pepsi_path + filename + " " + Data_path + flags + extra_flags
    os.system(command)
    with open('find_r0_dro.log', 'r') as f:
        lines = f.readlines()
    for line in lines:
        if "Best d_rho found" in line:
            line_tmp=line.split(':')[1]
            default_dro=float(line_tmp.split('e/A')[0])
        if "Best r0 found" in line:
            line_tmp=line.split(':')[1]
            default_r0=float(line_tmp.split('A')[0])
    
    # generate file with all fitted parameters
    params_file = 'params_' + str(j)
    with open(params_file,'w') as f:
        f.write('r0             \tdrho          \tI0          \tcst        \tchi2r\n')
    drho_list = []
    r0_list   = []
    cst_list  = []
    I0_list   = []
    for i in range(first_frame,last_frame+1):
        filename = "Free_AA_frame%dchain_%d.log" % (i,j)
        try:
            with open(filename, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if "Best d_rho found" in line:
                        line_tmp=line.split(':')[1]
                        value=float(line_tmp.split('e/A')[0])/default_dro
                        drho_list.append(value)
                        drho = str(value)
                    if "Best r0 found" in line:
                        line_tmp=line.split(':')[1]
                        value=float(line_tmp.split('A')[0])/default_r0
                        r0_list.append(value)
                        r0 = str(value)
                    if "Constant value" in line:
                        line_tmp=line.split(':')[1]
                        value=float(line_tmp)
                        cst_list.append(value)
                        cst = str(value)
                    if "I(0)" in line:
                        line_tmp=line.split(':')[1]
                        value=float(line_tmp)
                        I0_list.append(value)
                        I0 = str(value)
                    if "Chi^2" in line:
                        line_tmp=line.split(':')[1]
                        value=float(line_tmp)
                        chi2r=str(value)
        except:
            print(filename, 'parameter not found in log file - skipping it.')
        with open(params_file,'a') as f:
            f.write('%s \t%s \t%s \t%s \t%s\n' % (r0,drho,I0,cst,chi2r))

    # calculate mean parameter values and convert to string
    # for I0 the mean is the log-mean
    r0 = str(np.mean(r0_list))
    drho = str(np.mean(drho_list))
    logI = np.log10(I0_list)
    logI_mean = np.mean(logI)
    I0_mean = 10**logI_mean
    I0 = str(I0_mean)
    cst_mean = np.mean(cst_list)
    cst = str(cst_mean)

    # calculate chi2r for mean values
    chi2r_mean,Ifit_mean,q,Iexp,dI = get_chi2r(0,Pepsi_path,Data_path,flags,Ifit_column,first_frame,last_frame,r0,drho,I0,cst)
    with open(fit_file,'a') as f:
        f.write('#  \n')
        f.write('#  Fitting with mean parameter values:\n')
        f.write('#  I0_mean         \tcst                    \tchi2r \n')
        f.write('%s       \t%s  \t%1.2f    \n' % (I0,cst,chi2r_mean))

    #### fit parameters I0 and cst
    def func(q,S,B):
        I = S * Ifit_mean + B
        return I
    S0 = 1.00
    B0 = 0.00
    popt,pcov = curve_fit(func,q,Iexp,sigma=dI,p0=[S0,B0])
    S_opt = popt[0]
    B_opt = popt[1]

    I0_opt = str(I0_mean*S_opt)
    cst_opt = str(cst_mean+B_opt)
    
    with open(fit_file,'a') as f:
        f.write('#  \n')
        f.write('#  Fitting with refined parameter values:\n')
        f.write('#  I0               \tcst                \tchi2r \n')
     
    # calculate Ifit_opt and chi2r_opt
    Ifit_opt = func(q,*popt)
    R_opt = (Iexp - Ifit_opt)/dI
    nu = len(q)-4 # degrees of freedom with 4 being the number of fitted parameters
    chi2r_opt = sum(R_opt**2)/nu

    # write fitted intensity to file
    Ifit_file = 'Ifit_opt_%d.dat' % j
    with open(Ifit_file,'w') as f:
        f.write('#  I0_opt = %s\n' % I0_opt)
        f.write('#  chi2r_opt = %f\n' % chi2r_opt)
        f.write('#  \n')
        f.write('#  q            Iexp            dI              Ifit_mean       Ifit_opt\n')
        for i in range(len(q)):
            f.write('%e  \t%e  \t%e  \t%e  \t%e \n' % (q[i],Iexp[i],dI[i],Ifit_mean[i],Ifit_opt[i]))
    
    # write optimal parameters to fit file
    with open(fit_file,'a') as f:
        f.write('%s      \t%s    \t%1.2f    \n' % (I0_opt,cst_opt,chi2r_opt))

    # write chi2r and parameters to summary file
    with open(summary_file,'a') as f:
        f.write('%d          \t%f    \t%f    \t%s      \t%s   \t%s   \t%s    \t%f    \t%e\n' % (j,chi2r_mean,chi2r_opt,r0,drho,I0_opt,cst_opt,S_opt,B_opt))

print('\n\n Done with Pepsi part! \n')

# skip last (for SAXS data)

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
    outfile_exp = 'exp_data_' + str(combination) + '.dat'
    outfile_calc = 'calc_data_' + str(combination) + '.dat'

    # print header to exp file
    with open(outfile_exp, 'w') as f:
        f.write("# DATA=SAXS PRIOR=GAUSS\n")

    # make first column of header to calc file
    header = "# label"

    for j in range(len(D)):
  
        #make EXP file
        footerlines = 0
        if D[j] == -1: 
            data_filename = '/lindorffgrp-isilon/thomasen/hnRNPA1/FLA1_SAXSdata_08042020_averaged_newest/' + str(sys.argv[1]) + '.dat'
            headerlines = 0
        else:
            data_filename = ''
            headerlines = 1
        footerlines = 0
        q,Iexp,dI = np.genfromtxt(data_filename,skip_header=headerlines,skip_footer=footerlines,usecols=[0,1,2], unpack=True)
        if D[j] == -1:
            q = q
            Iexp = Iexp
            dI = dI
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
        frame_number = i+1 # index with 1 instead of 0 for frame numbers in calc file
        frame_line = "frame_%d" % frame_number
        for j in range(len(D)):
            headerlines = 6
            footerlines = 0
            Ifit_column = 4
            if D[j] == -1:
                contrast = 0
                Ifit_column = 3
            elif D[j] == 0:
                contrast = 1
            elif D[j] == 42:
                contrast = 2
            elif D[j] == 70:
                contrast = 3
            fit_filename = "Free_AA_frame%dchain_%d.fit" % (i,contrast)
            q,Iexp,dIexp,Ifit_free = np.genfromtxt(fit_filename,skip_header=headerlines,skip_footer=footerlines,usecols=[0,1,2,Ifit_column],unpack=True)
            for I_line in Ifit_free:
                frame_line += " \t %e" % I_line
        frame_line += '\n'  
        with open(outfile_calc,'a') as f:             
            f.write(frame_line)
    
print('\n Done with make_calc_exp part! \n')

#   This is a small script for reweighting molecular simulations using
#   the Bayesian/MaxEnt (BME) approach
#   Copyright (C) 2018 Sandro Bottaro (name . surname @ iit . it)
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License V3 as published by
#   the Free Software Foundation, 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function

import numpy as np
from scipy import optimize #version >0.13.0 for dogleg optimizer
import sys
import warnings
        
# these are known data types and error priors
exp_types = ["NOE","JCOUPLINGS","CS","SAXS","DIST","RDC"]
prior_types = ["GAUSS"]

# calculate relative entropy
def srel(w0,w1):
    idxs = np.where(w1>1.0E-50)
    return np.sum(w1[idxs]*np.log(w1[idxs]/w0[idxs]))

# back calculate from structure. Simple weighted average
def back_calc(sim_data,weights):
    return np.sum(weights[:,np.newaxis]*sim_data,axis=0)

# calculate chi square 
def chi_square(exp_data,sim_data,weights,bounds=None):
    diff = back_calc(sim_data,weights)-exp_data[:,0]
    if(bounds==None):
        return np.sum((diff*diff)/exp_data[:,1])
    else:
        idxs = [0.0 if(( diff[i] >= 0 and bounds[i][1] == 0) or (diff[i] <= 0 and bounds[i][0] == 0)) else 1.0 for i in range(len(bounds))]
        diff *= idxs                
        return np.sum((diff*diff)/exp_data[:,1])

# read experimental data
def read_exp(filename):
    
    fh = open(filename)
    first = fh.readline().split()
    
    # do some checks on first line
    assert first[0] == "#", "Error. First line of exp file %s must be in the format \# DATA=XXX PRIOR=XXX" % filename

    # find data type and data prior
    data_type = (first[1].split("=")[1]).strip()
    prior_type = (first[2].split("=")[1]).strip()
    assert data_type in exp_types , "Error. DATA must be one of the following: %s " % (exp_types)
    assert prior_type in prior_types , "Error. PRIOR must be one of the following: %s " % (prior_types)
    
    # noe power is 6 by default, but it can be changed in file 
    noe_power=6.0
    if(len(first)==4 and first[3]=="POWER"):
        noe_power = float((first[3].split("=")[1]).strip())
        assert noe_power>0., "# Error. POWER has to be larger than 0"
            
    ln = 0
    labels = []
    bounds = []
    exp_data = []
    
    for line in fh:
        # skip commented out data
        if("#" in line): continue

        ln += 1
        lsplit = line.split()
        if(len(lsplit) != 3 and len(lsplit) != 4):
            print("# ERROR: experimental line in %s not valid" % filename)
            print(line),
            sys.exit(1)
                
        labels.append( "%s"  % (lsplit[0]))
        
        # if lenght is 4, third and fourth columns are  r_avg and sigma
        avg0 = float(lsplit[1])
        sigma0 = float(lsplit[2])
        
        if(len(lsplit)==3):
            bounds.append([None,None])

        # if lenght is 4
        elif(len(lsplit)==4):

            # upper/lower bound
            if(lsplit[3] == "UPPER"):
                # upper/lower is inverted because of 1/r^n dependence of NOE
                if(data_type=="NOE"):
                    bounds.append([None, 0.0])
                else:
                    bounds.append([0.0,None])
            elif(lsplit[3] == "LOWER"):
                if(data_type=="NOE"):
                    bounds.append([0.0,None])
                else:
                    bounds.append([None, 0.0])


        if(data_type=="NOE"):
            avg = np.power(avg0,-noe_power)
            sigma = (noe_power*avg*sigma0/(avg0))**2
            exp_data.append([avg,sigma])
        else:
            exp_data.append([avg0,sigma0*sigma0])
                
    assert(len(labels) == len(bounds))
    assert(len(labels) == len(exp_data))
    
    return labels, bounds, exp_data, data_type, noe_power


def read_sim(filename):
    
    data = np.array([[float(x) for x in line.split()[1:]] for line in open(filename) if("#" not in line)])
    return data


def weight_exp(exp_file,sim_file,w0,w_opt,outfile,rows=[],cols=[]):
        
    fname_stats = "%s.stats.dat" % (outfile)
    fh_stats  = open(fname_stats,'w')

    labels, bounds, exp_data, data_type,noe_power =  read_exp(exp_file)
            
    sim_data = read_sim(sim_file)
    
    if(len(rows)==0): rows = range(sim_data.shape[0])
    if(len(cols)==0): cols = range(sim_data.shape[1])

    exp_data = np.array(exp_data)[cols,:]
    sim_data = np.array(sim_data)[rows,:]
    sim_data = sim_data[:,cols]
    bounds = [bounds[k] for k in cols]
    labels = [labels[k] for k in cols]
    assert len(w_opt) ==  len(sim_data),"# Number of weights %d not equal to number of data from simulation (%d) " % (len(w_opt),len(sim_data))

    if(data_type=="NOE"):
        sim_data = np.power(sim_data,-noe_power)

    L0 = 1.0
    L1 = 1.0
    if(data_type=="RDC"): 
        ee = np.array([e[0] for e in exp_data])
        prediction0 = np.average(sim_data,axis=0,weights=w0)
        prediction1 = np.average(sim_data,axis=0,weights=w_opt)
        L0 = np.dot(ee,prediction0)/np.dot(prediction0,prediction0)
        L1 = np.dot(ee,prediction1)/np.dot(prediction1,prediction1)
        print("# RDC scaling factors %8.4e, %8.4e" % (L0,L1))


    assert len(labels) == len(sim_data[0]), \
        "# Number of rows in %s not equal to number of columns in %s" % (len(labels),len(sim_data[0]))
    

        
    # before/after reweighting
    chisq_0 = chi_square(exp_data,L0*sim_data,w0,bounds=bounds)/len(exp_data)
    bcalc_0 = back_calc(L0*sim_data,w0)
        
    chisq_1 = chi_square(exp_data,L1*sim_data,w_opt,bounds=bounds)/len(exp_data)
    bcalc_1 = back_calc(L1*sim_data,w_opt)
        
    # write to file.
    fh_stats.write("# %s vs. %s srel=%8.4e\n" % (exp_file.split("/")[-1],sim_file.split("/")[-1],-srel(w0,w_opt)))
    fh_stats.write("#  %-15s %9s %9s %9s %9s \n" % ("Label","avg exp","sigma exp","avg before","avg after"))
    chisq_dist_0 = 0.0
    chisq_dist_1 = 0.0
    for l in range(exp_data.shape[0]):
        
        # convert to distances if is NOE
        fh_stats.write("   %-15s %9.3e %9.3e %9.3e %9.3e " % (labels[l],exp_data[l,0],np.sqrt(exp_data[l,1]),bcalc_0[l],bcalc_1[l]))
        if(data_type == "NOE"):
            exp_avg  = np.power(exp_data[l,0],-1./noe_power)
            sim_0  = np.power(bcalc_0[l],-1./noe_power)
            sim_1  = np.power(bcalc_1[l],-1./noe_power)
            sigma = np.sqrt(exp_data[l,1])*(exp_avg/(noe_power*exp_data[l,0]))
            diff_0 =  sim_0 - exp_avg
            diff_1 =  sim_1 - exp_avg
            ff_0 = 1.0
            ff_1 = 1.0
            if((diff_0 >= 0 and bounds[l][0] == 0) or (diff_0 <=0 and bounds[l][1] ==0 )):
                ff_0 = 0.0
            if((diff_1 >= 0 and bounds[l][0] ==0) or (diff_1 <=0 and bounds[l][1] ==0 )):
                ff_1 = 0.0
            chisq_dist_0 += ((diff_0)**2/(sigma*sigma))*ff_0
            chisq_dist_1 += ((diff_1)**2/(sigma*sigma))*ff_1
            fh_stats.write(" %9.3e %9.3e %9.3e %9.3e \n" % (exp_avg,sigma,sim_0,sim_1))
        else:
            fh_stats.write( "\n")

    fh_stats.write("#  %-5s %8.4f %8.4f \n" % ("chi2",chisq_0,chisq_1))
    if(data_type=="NOE"):
        chisq_dist_0 /= len(exp_data)
        chisq_dist_1 /= len(exp_data)
        fh_stats.write("#  %-5s %8.4f %8.4f \n" % ("chi2 distance space",chisq_dist_0,chisq_dist_1))
        fh_stats.close()
        return chisq_dist_0,chisq_dist_1
    else:
        fh_stats.close()
        return chisq_0,chisq_1
    

# reweight class. 
class Reweight:

    # initialize
    def __init__(self,verbose=True,w0=[],kbt=None, opt_method="L-BFGS-B"):

        self.exp_data = []
        self.sim_data = []
        self.labels = []
        self.bounds = []
        self.verbose = verbose
        
        self.result = None
        self.method = None
        self.lambdas = None
        self.success = None
        self.renormalize_w0 = False
        
        self.w0 = []
        self.w_opt = []
        if(len(w0)!=0):
            w0 = np.array(w0)
            if(kbt!=None):
                print("# Assuming weights given as minus free energies. w=exp(bias/kbt) kbt=%8.4f "  % kbt)
                w0 = np.exp((w0-np.max(w0))/kbt)
            print("# Set non-uniform initial weights from file. Sum=", np.sum(w0), len(w0))
            self.w0 = np.copy(w0)/np.sum(w0)
        self.opt_method = opt_method
    # add data to class
    def add_data(self,exp_data,sim_data,labels,bounds,cols):

        #if(len(rows)==0): rows = range(sim_data.shape[0])
        #if(len(cols)==0): cols = range(sim_data.shape[1])
        
        #sim_data1 = np.array(sim_data)[rows,:]
        #sim_data1 = sim_data1[:,cols]


        #print(exp_data)
        #print(np.copy(exp_data))
        if(len(self.sim_data)==0):
            self.sim_data = sim_data
            self.exp_data = np.copy(exp_data)[cols,:]
        else:
            self.sim_data = np.concatenate((self.sim_data,sim_data),axis=1)
            self.exp_data = np.concatenate((self.exp_data,np.array(exp_data)[cols,:]),axis=0)
        self.labels.extend([labels[k] for k in cols])
        self.bounds.extend([bounds[k] for k in cols])
        
    # load data from file
    def load(self,exp_file,sim_file,rows=[],cols=[]):

        # read experimental data 
        labels, bounds, exp_data, data_type,noe_power =  read_exp(exp_file)

        # read simulation data
        sim_data1 = read_sim(sim_file)
        
        if(len(rows)==0): rows = range(sim_data1.shape[0])
        if(len(cols)==0): cols = range(sim_data1.shape[1])
        if(len(self.w0)!=0 and self.renormalize_w0==False):
            self.w0 = np.array(self.w0)[rows]/np.sum(np.array(self.w0[rows]))
            self.renormalize_w0= True
            
        labels = np.array(labels)[cols]
        sim_data = np.array(sim_data1)[rows,:]
        sim_data = sim_data[:,cols]
        exp_data = np.array(exp_data)[cols, :]

        if(data_type=="NOE"):
            sim_data = np.power(sim_data,-noe_power)

        if(data_type=="RDC"):
            ee = np.array([e[0] for e in exp_data])
            if(len(self.w0)==0):
                prediction = np.average(sim_data,axis=0)
            else:
                prediction = np.average(sim_data,axis=0,weights=self.w0)
            L = np.dot(ee,prediction)/np.dot(prediction,prediction)
            print("# RDC scaling factor %8.4e" % L)
            sim_data *=	L
            
        assert len(labels) == sim_data.shape[1], "# Number of rows in exp (%d) not equal to number of columns in calc (%d)" % (len(labels),sim_data.shape[1])

        # now do some sanity checks on data 
        mins = np.min(sim_data,axis=0)
        maxs = np.max(sim_data,axis=0)
        norm_factor = []
        for l in range(len(exp_data)):

            if(bounds[l][0] == 0.0  and mins[l]> exp_data[l][0]):
                print("# Warning: expt boundary %20s=%-10.4e is smaller than minimum value in simulation %-10.4e"\
                      % (labels[l],exp_data[l][0],mins[l]))
            if(bounds[l][1] == 0.0  and maxs[l] < exp_data[l][0]):
                print("# Warning: expt boundary %20s=%-10.4e is larger than maximum value in simulation %-10.4e"\
                      % (labels[l],exp_data[l][0],mins[l]))
                
            if(exp_data[l][0] > maxs[l]):
                print("# Warning: expt average %20s=%-10.4e is larger than maximum value in simulation %-10.4e"\
                      % (labels[l],exp_data[l][0],maxs[l]))
            if(exp_data[l][0] < mins[l]):
                print("# Warning: expt average %20s=%-10.4e is smaller than minimum value in simulation %-10.4e"\
                      % (labels[l],exp_data[l][0],mins[l]))

            # if exp data is non-zero, normalize
            ff = 1.0
            
            if(np.abs(exp_data[l][0])>1.0E-05):
                ff = exp_data[l][0]
                exp_data[l][1] /= ff**2
                sim_data[:,l] /= ff
                exp_data[l][0] /= ff

        self.add_data(exp_data,sim_data,labels,bounds,cols)




    # Optimize
    def optimize(self,theta,method="MAXENT"):

        def func_maxent_gauss(lambdas):
            
            # weights
            arg = -np.sum(lambdas[np.newaxis,:]*self.sim_data,axis=1)
            # 
            overflow = np.zeros(arg.shape)
            tmax = np.log((sys.float_info.max)/50.)
            tmin = np.log((sys.float_info.min)/50.)
            overflow[np.where(arg>tmax)] =  tmax
            overflow[np.where(arg<tmin)] =  tmin
            arg -= overflow
            arg += overflow

            #if(np.max(arg)>300.):arg -= np.max(arg)
            ww  = (self.w0*np.exp(arg))/np.exp(overflow)
            #ww  = self.w0*np.exp(arg)
            # normalization 
            zz = np.sum(ww)
            if(zz < 1.0E-300):
                print("# ERROR. Sum of weights is very small (%4.2e)"  %(zz))
                print("# Perhaps you cannot fit the data. You may try with a larger theta?")
                print("# These are the top 10 frames")
                for j in np.argsort(ww[:10]):
                    print("%10d %10.4e" % (j,ww[j]))
                    
                sys.exit(1)
 
            ww /= zz
            # new averages
            avg = np.sum((ww[:,np.newaxis]*self.sim_data), axis=0)
            
            # errors are rescaled by factor theta
            err = (self.theta)*(self.exp_data[:,1])
            
            # gaussian integral
            eps2 = 0.5*np.sum((lambdas*lambdas)*err)
            # experimental value 
            sum1 = np.dot(lambdas,self.exp_data[:,0])
            fun = sum1 + eps2+ np.log(zz)
            
            # gradient
            jac = self.exp_data[:,0] + lambdas*err - avg
                   
            # divide by theta to avoid numerical problems
            return  fun/self.theta,jac/self.theta

        def hess_maxent_gauss(lambdas):
            
            # weights
            arg = -np.sum(lambdas[np.newaxis,:]*self.sim_data,axis=1)
            # 
            overflow = np.zeros(arg.shape)
            tmax = np.log((sys.float_info.max)/50.)
            tmin = np.log((sys.float_info.min)/50.)
            overflow[np.where(arg>tmax)] =  tmax
            overflow[np.where(arg<tmin)] =  tmin
            arg -= overflow
            arg += overflow

            #if(np.max(arg)>300.):arg -= np.max(arg)
            ww  = (self.w0*np.exp(arg))/np.exp(overflow)
            #ww  = self.w0*np.exp(arg)
            # normalization 
            zz = np.sum(ww)
 
            ww /= zz
            
            # errors are rescaled by factor theta
            err = (self.theta)*(self.exp_data[:,1])
                        
            # hessian
            q_w = np.dot(ww,self.sim_data)
            hess = np.einsum('k, ki, kj->ij',ww,self.sim_data,
                   self.sim_data) - \
                   np.outer(q_w,q_w) + \
                   np.diag(err)
                   
            # divide by theta to avoid numerical problems
            return  hess/self.theta

        
        
        def func_ber_gauss(w):
            
            bcalc = np.sum(w[:,np.newaxis]*self.sim_data,axis=0)
            diff = bcalc-self.exp_data[:,0]
            idxs = [0.0 if(( diff[i] >= 0 and self.bounds[i][1] == 0) or (diff[i] <= 0 and self.bounds[i][0] == 0)) else 1.0 for i in range(len(self.bounds))]
            diff *= idxs
            chi2_half =  0.5*np.sum(((diff**2)/(self.exp_data[:,1])))

            idxs = np.where(w>1.0E-50)
            srel = self.theta*np.sum(w[idxs]*np.log(w[idxs]/self.w0[idxs]))
            return chi2_half+srel

            
        self.theta = theta
        nf = self.sim_data.shape[0]
        # if weight are not set, use uniform
        if(len(self.w0)==0): self.w0 = np.ones(nf)/(1.*nf)
        self.method = method
        self.exp_data = np.array(self.exp_data)
        self.sim_data = np.array(self.sim_data)
        #print("# exp data: %s" % (str(self.exp_data.shape)))
        #print("# calc data: %s" % (str(self.sim_data.shape)))
        print("# theta: %s" % str(self.theta))
        # first, calculate initial chi squared and RMSD

        chi_sq0 = chi_square(self.exp_data,self.sim_data,self.w0,self.bounds)/len(self.exp_data)

        #print("Initial average chi square %10.4f, srel %10.4f " % (chi_sq0, srel(self.w0,self.w0)))
        
        if(method=="MAXENT"):
            
            opt={'maxiter':50000,'disp':False}
            meth = self.opt_method
            lambdas=np.zeros(self.exp_data.shape[0])

            result = optimize.minimize(func_maxent_gauss,lambdas,options=opt,method=meth,jac=True,bounds=self.bounds)
            #else:
             #   if meth == 'trust':
             #       meth = 'trust-exact'
             #   result = optimize.minimize(func_maxent_gauss,lambdas,
	#		        options=opt,method=meth, jac=True, 
#			        hess=hess_maxent_gauss)
            arg = -np.sum(result.x[np.newaxis,:]*self.sim_data,axis=1)
            if(np.max(arg)>300.):arg -= np.max(arg)

            w_opt = self.w0*np.exp(arg)
            w_opt /= np.sum(w_opt)
            self.lambdas = np.copy(result.x)
            
        if(method=="BER"):
            opt={'maxiter':2000,'disp': self.verbose,'ftol':1.0e-20}
            cons = {'type': 'eq', 'fun':lambda x: np.sum(x)-1.0}
            bounds = [(0.,1.)]*len(self.w0)  
            meth = "SLSQP"
            print("# Bayesian Ensemble Refinement. Useful for testing purposes and ")
            print("# when the number of experimental data is larger than the number of samples.")
            
            result = optimize.minimize(func_ber_gauss,self.w0,constraints=cons,options=opt,method=meth,bounds=bounds)

            w_opt = result.x

        self.success = result.success
        
        if(result.success):
            chi_sq1 = chi_square(self.exp_data,self.sim_data,w_opt,self.bounds)/len(self.exp_data)
            self.w_opt = np.copy(w_opt)
            srel1 = srel(self.w0,w_opt)
            #print("# Minimization successful")
            #print("Final average chi squared   %10.4f, srel %10.4f " % (chi_sq1, srel(self.w0,w_opt)))
            return chi_sq0, chi_sq1, -srel1 
        else:
            print("# ERROR - Minimization failed. Perhaps theta is too small, or there is some bug.")
            print("# exit message:", result.message)
            sys.exit(1)
        return 0
    
    def get_lambdas(self):

        assert self.method=="MAXENT", "# Error. You can get lambdas only if method=MAXENT"
        assert self.lambdas!= None, "# Error. Minimization not successful, no lagrange multipliers"
        return self.lambdas
    
    def get_weights(self):
        
        assert self.success == True
        return self.w_opt

        
    def weight_exp(self,exp_file,sim_file,outfile,rows=[],cols=[]):
        
        return weight_exp(exp_file,sim_file,self.w0,self.w_opt,outfile,rows,cols)
        





'''        
        def func_maxent_laplace(lambdas):

            # weights
            ww  = self.w0*np.exp(-np.sum(lambdas[np.newaxis,:]*self.sim_data,axis=1))
            # normalization 
            zz = np.sum(ww)
            ww /= zz
            # new averages
            avg = np.sum((ww[:,np.newaxis]*self.sim_data), axis=0)
            # errors are rescaled by factor theta
            err = theta*self.exp_data[:,1]
            
            # integral error
            eps2 = np.sqrt(2.*err)/(1.-0.5*(err*lambdas**2))
            eps2 = np.sum(np.log(eps2))
    
            # experimental value 
            sum1 = np.dot(lambdas,self.exp_data[:,0])
            fun = sum1 + eps2+ np.log(zz)
            # gradient
            lap0 = -lambdas*err
            lap1 = lap0/(1.+0.5*(lambdas*lap0))
            jac = exp_avg[:,0] - lap1 - avg

            return  fun,jac
'''
            # write weights to file
            #fname_weights = "%s.weights.dat" % (self.outfile)
            #fh_weights  = open(fname_weights,'w')
            #fh_weights.write("# Reweighting. n_data=%d, n_samples=%d \n" % (self.exp_data.shape[0],self.sim_data.shape[0]))
            #fh_weights.write("# %-5s %8.4f  \n" % ("Theta",self.theta))
            #neff1 = np.exp(np.sum(-self.w_opt*np.log(self.w_opt/self.w0)))
            #fh_weights.write("# %-5s %8.4f  \n" % ("neff",neff1))
            #fh_weights.write("# %-5s %8.4e %8.4e \n" % ("sum weights",np.sum(self.w0),np.sum(self.w_opt)))
            #fh_weights.write("# %8s %8s \n" % ("before","after"))
            #for l in range(len(self.w0)):
            #    fh_weights.write("  %10.5e %10.5e \n" % (self.w0[l],self.w_opt[l]))
            #fh_weights.close()

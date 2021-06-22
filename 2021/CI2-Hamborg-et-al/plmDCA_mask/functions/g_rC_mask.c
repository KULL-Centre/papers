#include <math.h>
#include "mex.h"

/*
This is a reworked version of the pseudolikelihood objective contained in Mark Schmidt's thesis code (http://www.di.ens.fr/~mschmidt/Software/thesis.html). The copyright conditions for the original code is included below.
---------------------------------

Copyright 2005-2012 Mark Schmidt. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ''AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* variables */
    char param;
    int i, s, t, n, nInstances, nNodes, nStates,
            *y, y1, y2, *rint,r,*mask;
    double *weights, *grad1, *grad2, *logPot, *z, *fval, *h_r, *J_r, *nodeBel, *lambdas;
    
    /* input */
    y = (int*)mxGetPr(prhs[0]);
    weights = mxGetPr(prhs[1]);
    h_r = mxGetPr(prhs[2]);
    J_r = mxGetPr(prhs[3]);
    lambdas = mxGetPr(prhs[4]);
    rint = (int*)mxGetPr(prhs[5]);
    mask = (int*)mxGetPr(prhs[6]);
    
    /* compute sizes */
    nInstances = mxGetDimensions(prhs[0])[0];
    nNodes = mxGetDimensions(prhs[0])[1];
    nStates = mxGetDimensions(prhs[2])[1];
    
    /* allocate memory */
    logPot = mxCalloc(nStates, sizeof(double));
    z = mxCalloc(1, sizeof(double));
    nodeBel = mxCalloc(nStates, sizeof(double));
    
    /* output */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    fval = mxGetPr(plhs[0]);
    *fval = 0;
    plhs[1] = mxCreateDoubleMatrix(nStates, 1, mxREAL);
    grad1 = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(nStates*nStates*(nNodes-1), 1, mxREAL);
    grad2 = mxGetPr(plhs[2]);
    
    for (i = 0; i < nStates*nStates*(nNodes-1); ++i)
    {   
        if(grad2[i]!=0){
            printf("grad2 = %f\n",grad2[i]);
        }

    }

    r=rint[0]-1;


    for(i=0;i < nInstances;i++) {
	/*Some notes on variable names:
	logPot contains, for the current sequence i, the exponentials in the pseudolikelihood: logPot(s)=h_r(s)+sum_{j!=r}J_{rj}(s,sigma^(i)_j).
	nodeBel is the conditional probability P(sigma_r=s|sigma_{\r}=sigma^(i)_{\r}), i.e., nodeBel(s) = e^[ logPot(s) ] / sum_l e^[ logPot(l) ].
	z is the denominator of nodeBel.*/

        for(s=0;s < nStates;s++) {
            logPot[s] = h_r[s];
        }

        for(n = 0;n < nNodes;n++) {
            if(n!=r && mask[n]!=0) {
	        y2 = y[i + nInstances*n];       
                for(s=0; s<nStates; s++) {
                    logPot[s] += J_r[s+nStates*(y2+nStates*(n-(n>r)))];               
                }
            }
            
        }
      
        z[0] = 0;
        for(s = 0; s < nStates; s++) {
            z[0] += exp(logPot[s]);
        }
        *fval -= weights[i]*logPot[y[i+nInstances*r]];
        *fval += weights[i]*log(z[0]);
 
        
        
	/*Gradient:*/
        
 
        for(s = 0; s < nStates; s++) {
            nodeBel[s] = exp(logPot[s] - log(z[0]));
        }
                       
        y1 = y[i + nInstances*r]; 
        grad1[y1] -= weights[i]*1;
            
        for(s=0; s < nStates; s++) {
            grad1[s] += weights[i]*nodeBel[s];
        }
   
        for(n=0;n<nNodes;n++) {
            if(n!=r && mask[n]!=0)
            /*if(n!=r)*/ {
                y2 = y[i + nInstances*n];

                grad2[y1+nStates*(y2+nStates*(n-(n>r)))] -= weights[i];                

                for(s=0;s<nStates;s++) {
                   grad2[s+nStates*(y2+nStates*(n-(n>r)))] += weights[i]*nodeBel[s];
                }	
            }  
        }
        
    }

  
    
    
    
    
    /*Add contributions from R_l2*/

    for(s = 0; s < nStates; s++) {
        *fval += lambdas[0]*h_r[s]*h_r[s];
        grad1[s] += lambdas[0]*2*h_r[s]; 
    }

    for(n = 0;n < nNodes;n++) {
        if(n!=r) {
            for(s = 0; s < nStates; s++) {
                for(t = 0; t < nStates; t++) {                   
                        *fval += lambdas[1]*J_r[s+nStates*(t+nStates*(n-(n>r)))]*J_r[s+nStates*(t+nStates*(n-(n>r)))];
                        grad2[s+nStates*(t+nStates*(n-(n>r)))] += lambdas[1]*2*J_r[s+nStates*(t+nStates*(n-(n>r)))];                    
                }
            }
        }
    }
             
    mxFree(logPot);
    mxFree(z);
    mxFree(nodeBel);
    return;
}

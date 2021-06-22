UPDATES:
v2:
-Regularization strengths are now chosen automatically by the program, i.e. they are no longer input arguments. If the number of nonredundant homogous sequences B_eff is >500, the previously recommended standard values lambda_h=lambda_J=0.01 are used, whereas for B_eff<=500 the value is now instead taken as 0.1-(0.1-0.01)*B_eff/500. This because we have observed a slight boost in accuracy when using stronger regularization at low B_eff.
-Contributions from the gap-state are now excluded during computation of the coevolution score. This has been found to yield a small but consistent improvement.

---------------------------------------------
This is a MATLAB-implementation of the asymmetric version of plmDCA.
The program uses 'minFunc' by Mark Schmidt (contained in the folder '3rd_party_code'), which can be found at http://www.di.ens.fr/~mschmidt/Software/minFunc.html.


HOW TO USE:
	1) In MATLAB, go to the directory containing mexAll.m, and give the command 
	>mexAll 
	to create mex files appropriate for your system.
	2) Call the program as
	plmDCA_asymmetric(fastafile,outputfile,reweighting_threshold,nr_of_cores)

 
INPUTS TO PLMDCA: 
	fastafile:
	Alignment file in FASTA format. Inserts, which should be represented by '.' and lower-case letters (as is standard in the Pfam download), are removed automatically by the program.

	outputfile:
	This becomes a file with N(N-1)/2 rows (N=domain length), each with three entries: residue i, residue j, and interaction score for the pair (i,j).

	reweighting_threshold:
	Required fraction of nonidentical amino acids for two sequences to be counted as independent (typical values: 0.1-0.3).
	Note that this is not the threshold 'x' as defined in the papers, but 1-x.
 	
	nr_of_cores:
	The number of cores to use on the local machine. 
	If this argument is >1, the program calls functions from MATLAB's Parallel Computing Toolbox.


TYPICAL CALL (ON A QUAD-CORE MACHINE):
	plmDCA_asymmetric('PF00014_full.txt','PF00014_scores.txt',0.2,4)



---------------------------------------------
 Copyright 2014 - by Magnus Ekeberg (magnus.ekeberg@gmail.com)
 All rights reserved
 
 Permission is granted for anyone to copy, use, or modify this
 software for any uncommercial purposes, provided this copyright 
 notice is retained, and note is made of any changes that have 
 been made. This software is distributed without any warranty, 
 express or implied. In no event shall the author or contributors be 
 liable for any damage arising out of the use of this software.
 
 The publication of research using this software, modified or not, must include 
 appropriate citations to:

 	M. Ekeberg, C. Lövkvist, Y. Lan, M. Weigt, E. Aurell, Improved contact
 	prediction in proteins: Using pseudolikelihoods to infer Potts models, Phys. Rev. E 87, 012707 (2013)

	M. Ekeberg, T. Hartonen, E. Aurell, Fast pseudolikelihood
	maximization for direct-coupling analysis of protein structure
	from many homologous amino-acid sequences, J. Comput. Phys. 276, 341-356 (2014)

# Data files related to BME reweighting

For thorough introduction and protocol for BME reweighting see (1).


## Files:
    - simulation_HN_NOE.dat:
        contains the backcalculated HN-NOEs for each 1 ns frame in the simulation. Rows = timestep/frame (ns),  columns = atompair in HN-NOE. First column is timestep/frame. Columns are orderede as in experimental file.  
    - simulation_methyl_NOE.dat
        same as above but for methyl NOE pairs
    - simulation_SAXS.dat: 
        contains the backcalculated scatter (SAXS) intensity for each 1 ns frame in the simulation. Rows = timestep/frame (ns), columns = q_i 
    - BME_simulation_weights.dat
        contains the optimized weights for each simulation frame from the BME method. 
    - BME_stats_SAXS.dat
        SAXS from simulation  after BME reweight
    - BME_stats_methyl_NOE.dat
         Methyl NOE from simulation  after BME reweight
     - hyperparameter_theta.dat
        contains calculated chi^2 vs Srel for determination of the hyperparameter theta in the BME method. 
    - exp_HN2_NOE.dat
        HN NOE experimental file from (2), excluding any intraresidues. 
    - exp_methyl_NOE.dat
        NOE experimental file from (2), excluding any intraresidues.
    - the experimental SAXS file used for reweighting can be found in SAXS_and_SANS directory. 
        


## CITATIONS
(1)
@article{BME,
author = {Bottaro S., Bengtsen T., Lindorff-Larsen K},
title = {Integrating Molecular Simulation and Experimental Data: A Bayesian/Maximum Entropy Reweighting Approach}, 
booktitle = {Structural Bioinformatics. Methods in Molecular Biology},
year = {2020}}

(2) 
@article{NMR_Nanodisc,
title = {{Solution structure of discoidal high-density lipoprotein particles with a shortened apolipoprotein A-I}},
author = {Bibow, Stefan and Polyhach, Yevhen and Eichmann, C{\'{e}}dric and Chi, Celestine N and Kowal, Julia and Albiez, Stefan and McLeod, Robert A and Stahlberg, Henning and Jeschke, Gunnar and G{\"{u}}ntert, Peter and Riek, Roland},
journal = {Nature Structural {\&} Molecular Biology},
year = {2017}}



## Supplementary Information scripts for Cagiada et al(2020) ##

This repository contains all the necessary notebooks and scripts to run the analysis for both of proteins described in the manuscript.
The notebooks should produce the figures that are presented in the manuscript.  
The static images of the proteins in the manuscript rendered can be visualized interactively using Chimera (https://www.cgl.ucsf.edu/chimera/) and the relevant scripts contained in chimera_script folder.

## NOTEBOOK INSTRUCTIONS ##

1. Download the jupyter notebooks in this repository.

2. Download the Input_dataset folder (which contain all the computational data needed to run the notebooks) and place it in the same folder as the notebooks.
   
3. Download the data from the four MAVEs from SI of the papers used as data sources in this work:
  - Suiter, et al, 2020 (https://doi.org/10.1073/pnas.1915680117)
  - Mighell, et al, 2018 (https://doi.org/10.1016/j.ajhg.2018.03.018)
  - Matreyek, et al, 2018 (https://doi.org/10.1038/s41588-018-0122-z)
   and place them in the Input_dataset folder (NB. If errors occur, please check if the filenames are the same as in the notebooks)
   
4. Create a folder called Figures in the same folder as the notebook.

5. Now you can open and run the jupyter notebooks.

## CHIMERA scripts ##

The chimera_scripts folder contains the chimera python scripts  necessary to create the protein renderings shown in the paper.
Before running them you should download the structures (PDB entry 5LPG for NUDT15 and 1D5R for PTEN) from the PDB.

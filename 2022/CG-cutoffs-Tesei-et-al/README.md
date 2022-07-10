[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6815068.svg)](https://doi.org/10.5281/zenodo.6815068)

# Improved Predictions of Phase Behaviour of Intrinsically Disordered Proteins by Tuning the Interaction Range

This repository contains Python code, [Jupyter](http://jupyter.org) Notebooks and simulation data for reproducing the work of the scientific paper _Improved Predictions of Phase Behaviour of Intrinsically Disordered Proteins by Tuning the Interaction Range_ by G. Tesei and K. Lindorff-Larsen.

### Layout

- `analyses.ipynb` Jupyter Notebook to analyze all the simulation data and generate plots
- `calc_conc.ipynb` Jupyter Notebook to calculate _c<sub>sat</sub>_ and _c<sub>con</sub>_ from direct-coexistence molecular simulations
- `prior.ipynb` Jupyter Notebook to carry out the analysis of the hydrophobicity scales collected by Simm et al. [DOI: 10.1186/s40659-016-0092-5](https://doi.org/10.1186/s40659-016-0092-5)
- `optimization/` Data and Python code related to the optimization of the residue-specific ``stickiness'' parameters 
- `SC/` Data and Python code related to single-chain simulations of the CALVADOS model. Simulations are performed using [HOOMD-blue](https://hoomd-blue.readthedocs.io/en/latest/) v2.9.3 installed with the [mphowardlab/azplugins](https://github.com/mphowardlab/azplugins)
- `MC/` Data and Python code related to multi-chain simulations of the CALVADOS model in slab geometry. Simulations are performed using [openMM](https://openmm.org/) v7.5

Further usage examples of the CALVADOS model are available at [KULL-Centre/CALVADOS](https://github.com/KULL-Centre/CALVADOS).

### Usage

To open the Notebooks, install [Miniconda](https://conda.io/miniconda.html) and make sure all required packages are installed by issuing the following terminal commands

```bash
    conda env create -f environment.yml
    source activate calvados
    jupyter-notebook
```

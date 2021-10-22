[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5005953.svg)](https://doi.org/10.5281/zenodo.5005953)

# CG model of liquid-liquid phase behaviour of IDPs

This repository contains Python code, [Jupyter](http://jupyter.org) Notebooks and simulation data for reproducing the work of the scientific paper _Accurate model of liquid-liquid phase behaviour of intrinsically-disordered proteins from optimization of single-chain properties_ by G. Tesei, T. K. Schulze, R. Crehuet, and K. Lindorff-Larsen [DOI: 10.1101/2021.06.23.449550](https://doi.org/10.1101/2021.06.23.449550).

### Layout

- `analyses_and_figures.ipynb` Jupyter Notebook to analyze all the simulation data and generate plots
- `calcKd.ipynb` Jupyter Notebook to calculate _B<sub>22</sub>_ and _K<sub>d</sub>_ from two-chain simulations
- `calcConc.ipynb` Jupyter Notebook to calculate _c<sub>sat</sub>_ and _c<sub>con</sub>_ from multi-chain molecular simulations in slab geometry
- `analysis_HP_scales.ipynb` Jupyter Notebook to carry out the analysis of the hydrophobicity scales collected by Simm et al. [DOI: 10.1186/s40659-016-0092-5](https://doi.org/10.1186/s40659-016-0092-5)
- `single-chain/code/` Python code to simulate and analyze simulations of a single IDP of a given sequence using [HOOMD-blue](https://hoomd-blue.readthedocs.io/en/latest/) 
- `two-chain/code/` Python code to perform two-chain simulations and trajectory analyses of the optimized CG-IDPs model using HOOMD-blue
- `optimization/code/` Python code and bash script to optimize the CG-IDPs model against experimental radii of gyration and intramolecular NMR PRE data
- `multi-chain/code/` Python code to simulate (using HOOMD-blue or openMM) and analyze multi-chain simulations of the CG-IDPs model in slab geometry

### Usage

To open the Notebooks, install [Miniconda](https://conda.io/miniconda.html) and make sure all required packages are installed by issuing the following terminal commands

```bash
    conda env create -f environment.yml
    source activate cg-idps
    jupyter-notebook
```

This repository also contains an unmodified copy of [PULCHRA PowerfUL CHain Restoration Algorithm Version 3.06](https://www.pirx.com/pulchra/) (archive file `pulchra_306.tgz`). 

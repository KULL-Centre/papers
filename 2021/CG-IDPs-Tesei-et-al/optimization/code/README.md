### Layout

- `simulate.py` Python code to run a Langevin dynamics simulation of an IDP at infinite dilution using [HOOMD-blue](https://hoomd-blue.readthedocs.io/en/latest/)
- `pulchra.py` Python code to covert CG trajectories into all-atom trajectories using [PULCHRA](https://doi.org/10.1002/jcc.20906) 
- `analysis.py` helper functions to set up simulations, analyse trajectories and generate pandas DataFrames containing protein sequences and experimental conditions
- `optimize.py` Python code to perform the Bayesian parameter-learning procedure and optimize the hydropathy parameters against experimental gyration radii and NMR PRE data
- `optimize.sh` SLURM job script to run the complete pipeline of the optimization procedure

### Usage

To run the optimization procedure, install [Miniconda](https://conda.io/miniconda.html) and make sure all the required packages are installed by issuing the following terminal commands

```bash
    conda env create -f environment.yml
    source activate cg-idps
```

where `environment.yml` can be found in the parent directory of this repository.

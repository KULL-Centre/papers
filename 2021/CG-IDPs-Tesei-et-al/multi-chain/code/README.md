### Layout

- `simulate.py` Python code to run a Langevin dynamics multi-chain simulation using [HOOMD-blue](https://hoomd-blue.readthedocs.io/en/latest/)
- `submit.py` Python script to generate and submit the SLURM job script to run multi-chain simulations
- `maps_calc.py` Python code to calculate non-electrostatic interaction energies and energy maps
- `calc_submit.py` Python script to generate and submit the SLURM job script to run maps_calc.py
- `analysis.py` helper functions including initProteins() which generates a pandas DataFrame containing sequences and solution conditions

### Usage

To run a single-chain simulation, install [Miniconda](https://conda.io/miniconda.html) and make sure all the required packages are installed by issuing the following terminal commands

```bash
    conda env create -f environment.yml
    source activate cg-idps
```

where `environment.yml` can be found in the parent directory of this repository.


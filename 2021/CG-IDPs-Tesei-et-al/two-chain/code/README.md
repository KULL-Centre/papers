### Layout

- `simulate.py` Python code to run a Langevin dynamics simulation of two IDPs using [HOOMD-blue](https://hoomd-blue.readthedocs.io/en/latest/)
- `submit.py` Python script to generate and submit the SLURM job script to run two-chain simulations
- `maps_calc.py` Python code to calculate non-electrostatic interaction energies and energy maps
- `maps_submit.py` Python script to generate and submit the SLURM job script to run maps_calc.py
- `rdf_calc.py` Python code to calculate RDFs of the mass centers of the chains
- `rdf_submit.py` Python script to generate and submit the SLURM job script to run rdf_calc.py
- `pulchra.py` Python code to convert two-chain trajectories from CG to all-atom
- `pulchra_submit.py` Python script to generate and submit the SLURM job script to run pulchra.py
- `calcPREs.py` Python code to calculate NMR PRE data from two-chain all-atom trajectories
- `calcPREs_submit.py` Python script to generate and submit the SLURM job script to run calcPREs.py
- `optTauc.py` Python code to identify the _𝜏<sub>c</sub>_ value minimizing the discrepancy between experimental and calculated intermolecular NMR PRE rates
- `optTauc_submit.py` Python script to generate and submit the SLURM job script to run optTauc.py
- `analysis.py` helper functions including initProteins() which generates a pandas DataFrame containing sequences and solution conditions

### Usage

To run a single-chain simulation, install [Miniconda](https://conda.io/miniconda.html) and make sure all the required packages are installed by issuing the following terminal commands

```bash
    conda env create -f environment.yml
    source activate cg-idps
```

where `environment.yml` can be found in the parent directory of this repository.


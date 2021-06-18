### Layout

- `simulate.py` Python code to run a Langevin dynamics simulation of an IDP at infinite dilution using [HOOMD-blue](https://hoomd-blue.readthedocs.io/en/latest/)
- `submit.py` Python script to (i) generate a pandas DataFrame containing sequences and solution conditions, (ii) generate the SLURM job script, (iii) submit single-chain simulations.
- `calc.py` Python code to calculate the gyration radius, scaling exponent and non-electrostatic interaction energy maps
- `calc_submit.py` Python script to generate and submit the SLURM job script to run calc.py

### Usage

To run a single-chain simulation, install [Miniconda](https://conda.io/miniconda.html) and make sure all the required packages are installed by issuing the following terminal commands

```bash
    conda env create -f environment.yml
    source activate cg-idps
```

where `environment.yml` can be found in the parent directory of this repository.

Insert your FASTA sequence and solution conditions in initProteins(), select the sequences to simulate and submit the SLURM job script by issuing the following terminal command

```bash
    python submit.py
```

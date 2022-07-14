Scripts and data to repeat the analyses in Blaabjerg et al.: 
"Rapid protein stability prediction using deep learning representations"

## Code
Overview of files:<br>
`src/run_pred.py` contains a pipeline for making ddG ML predictions.<br/>
`src/run_pred_SLURM.sh` contains a recommended SLURM script.<br/>

## Installation
(Tested on Linux using Miniconda).

1. Clone the repository.

2. Install and activate conda environment with requirements:<br> 
`conda create --name cavity-model python=3.6`<br>
`conda activate cavity-model`<br>
`conda install pandas scipy numpy=1.17.3 scikit-learn mpl-scatter-density pdbfixer=1.5 pytorch=1.2.0 cudatoolkit=10.0 biopython=1.72 openmm=7.3.1 matplotlib=3.1.1 seaborn ptitprince -c omnia -c conda-forge -c anaconda -c defaults`

3. Install reduce in the right directory. This program is used by the parser to add missing hydrogens to the proteins.<br/>
`cd src/pdb_parser_scripts`<br/>
`git clone` https://github.com/rlabduke/reduce.git <br/>
`cd reduce/`<br/>
`make`; `make install` # This might give an error but provides the reduce executable in this directory<br/>

## Execution
1. Execute the pipeline on using the python script `src/run_pipeline.py`.

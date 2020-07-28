# Computing, analyzing and comparing the radius of gyration and hydrodynamic radius in conformational ensembles of intrinsically disordered proteins

This is a Python code to perform calculation of the radius of gyration and hydrodynamic radius from a structual ensemble. This code also include a how-to guide to perform BME reweighting and re-generate figures and data in main text. Here we use Ensemble 1 from the paper as an exemplar.  

### Getting Started

The folder `/notebook` contains a guide to generate data and figures from the main text in jupyter notebook format.

The folder `/data` contains the input files for the notebook.

The conformational ensembles that were used to generate the data are available at https://doi.org/10.5281/zenodo.3964622


### Requirements

1. Python Jupyter and Matplotlib with Python >=3.4  https://jupyter.org/
2. MDtraj  http://mdtraj.org/1.9.3/
3. Numpy libraries (native in Python Jupyter)
4. BME  https://github.com/KULL-Centre/BME
5. uncertainties 3.1.2  https://pypi.org/project/uncertainties/
6. Pepsi-SAXS  https://team.inria.fr/nano-d/software/pepsi-saxs



### References 
You may consider reading and citing the following relevant papers as well:
```
@article{bottaro2018integrating,
  title={Integrating Molecular Simulation and Experimental Data: A Bayesian/Maximum Entropy Reweighting Approach},
  author={Bottaro, Sandro and Bengtsen, Tone and Lindorff-Larsen, Kresten},
  journal={bioRxiv},
  pages={457952},
  year={2018},
  publisher={Cold Spring Harbor Laboratory}
}
```
```
@article{nygaard2017efficient,
  title={An efficient method for estimating the hydrodynamic radius of disordered protein conformations},
  author={Nygaard, Mads and Kragelund, Birthe B and Papaleo, Elena and Lindorff-Larsen, Kresten},
  journal={Biophysical journal},
  volume={113},
  number={3},
  pages={550--557},
  year={2017},
  publisher={Elsevier}
}
```
### Acknowledgments


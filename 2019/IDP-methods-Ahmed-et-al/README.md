# Analyzing and comparing the radius of gyration and hydrodynamic radius in conformational ensembles of intrinsically disordered proteins

This is a Python code to perform calclation of the radius of gyration and hydrodynamic radius from a structual ensemble. This code also include a how-to guide to perform BME reweighting and re-generate figures and data in main text. Here, We use Ensemble 1 from the paper as an exemplar.  

### Getting Started

The folder `/notebook` contains a guide to generate data and figures from the main text in jupyter notebook format.

The folder `/data` contains the input files for the notebook.


### Requirements

1. Python >=3.4
2. BME 
3. Numpy, Scipy libraries
4. Jupyter and Matplotlib
5. MDtraj

### References 
You may consider reading and citing the following relevant references as well:
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


# Discovering functionally important sites in proteins

This repository contains Python code, [Jupyter](http://jupyter.org) Notebooks and used data for reproducing the work of the scientific paper _Discovering functionally important sites in proteins_ by M. Cagiada, S. Bottaro, S. Lindemose, S. M. Schenstrøm, A. Stein, R. Hartmann-Petersen and K. Lindorff-Larsen, DOI: to add
### Layout

- `Functional_site_model.ipynb` Jupyter Notebook to process and reproduce all the prediction uses in the paper and make new predictions on novel data.
- `catboost_model` folder which contains the trained Catboost model.
- `scores_GEMME` folder which contains the GEMME scores predicted for the he proteins used in the manuscript.
- `scores_rosetta` folder which contains the rosetta stability predictions predicted scores for the proteins used in the manuscript.
- `pdbs` folder which contains pdbs used for the proteins used in the manuscript.

### Usage

To run the Notebooks, we suggest to open it in [Google Colaboratory]([https://colab.research.google.com/).
To run the Notebook is necessary to follow different steps:
1 Run the first two cells, which install all the dependencies and load the required functions to run the predictor.
2 Run in the order all the cells marked as 'Preliminary operations' which allow to upload all the necessary information for running the notebook:
  - adding the target protein sequence
  - adding the file with the thermodynamic stability measure (for the file format look in the following section).
  - adding the file with the evolutionary conservation measure made (for the file format look in the following section).
  - adding the target protein pdb (which should cover all the target sequence or part of it)
3 Run the cell 'Variant prediction and residue classification' to load the model and using it to make predictions with the provided input informations.
4 Run the cell 'Download predictions'  which allow to generate and create the files with prediction at variant and residue level.


### File format:
1 All the input data should follow the same numeration of the uploaded target sequence. In case of different numeration, for themodynamic stability changes and conservation measures is it possible to input the difference between the uploaded target sequence and the file's sequence. This option is not available in the pdb uploading cell.

2 Input files containing thermodynamic stability changes prediction should be properly formatted.

----------------------------------
E.g. for Protein with 45 residues:

Mutation  Score[in kcal/mol]
M1A       2.4
M1C       1.2
...
M1W       1.3
C2A       0.2
..
..
Y45W       0.3

----------------------------------
Synonymous mutations should be skipped or reported as mutation to =. Stop mutations has to been removed from the list.
Examples of predictions can be found in 'score_rosetta' folder.

3 Input files containing evoluationary conservation measures should follow the same format rules as shown for stability prediction or the GEMME website output folder (remember to select in the notebook the format used). All the prediction should be generated using GEMME.

4 PBD files should follow the same numeration and should not contain extra residue/regions compared to the uploaded target sequence. The target structure has to be named as chain A.

### Output:
Running  the 'Download predictions' cell generates two output files, which can be downloaded as archive or display online in the colaboratory.
- The first file (named: 'prefix_variant_features.txt') contain a summary of the input features for each variant of the target sequence.
- The second file (named: 'prefix_variant_prediction.txt') contains the model predictions for each single variant of the target protein.
- The third file (named: 'prefix_residue_prediction.txt') instead contains information about the predicted residue class and statistic about variant classes for each residue.


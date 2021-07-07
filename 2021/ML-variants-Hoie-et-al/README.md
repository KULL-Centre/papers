## Supplementary Information scripts for ML-variants Hoie et al (2021) ##

This repository contains all the necessary notebooks and scripts to run the analysis described in the manuscript.

## DOWNLOAD AND SETUP ENVIRONMENT

```bash
svn export https://github.com/KULL-Centre/papers/trunk/2021/ML-variants-Hoie-et-al
cd ML-variants-Hoie-et-al
unzip data.zip

conda create --name ML-variants --file environment.yml
conda activate ML-variants
```

## NOTEBOOKS

- figures.ipynb: generates the main manuscript figures
- supp.ipynb: generates supplementary figures.

## PRE-PROCESSING AND TRAINING

```bash
# Merge datasets in data/raw folder
python src/utilities/prism_folder_merge.py

# Create new pre-processed .pkl file for random-forest script
python src/utilities/generate_fresh_dataset_current.py

# Train new random-forest models with output predictions in data/runs
bash train.sh
```

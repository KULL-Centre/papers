#!/bin/bash

echo "Generating new merged prism files"
python src/utilities/prism_folder_merge.py

echo "Generating new pythonic dataset"
python src/utilities/generate_fresh_dataset_current.py

#echo "Running ML models on dataset"
#cat train.sh | parallel -j 4

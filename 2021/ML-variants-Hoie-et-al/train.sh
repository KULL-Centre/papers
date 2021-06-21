#!/bin/bash

# Exclude all missing gemme and ros
for i in $(seq 0 39); do python -u "src/RandomForest_model.py" -x 2 -d "data/preprocessed.pkl" -i [$i] -n "paper__current__exclude2__ranknorm__all" ; done

for i in $(seq 0 39); do python -u "src/RandomForest_model.py" -x 2 -d "data/preprocessed.pkl" -i [$i] -n "paper__current__exclude2__ranknorm__ddg_dde_only" -f '["ros_aa_wt_p$", "gemme_aa_wt_p$"]' ; done

for i in $(seq 0 39); do python -u "src/RandomForest_model.py" -x 2 -d "data/preprocessed.pkl" -i [$i] -n "paper__current__exclude2__ranknorm__ddg_only" -f '["ros_aa_wt_p$"]' ; done

for i in $(seq 0 39); do python -u "src/RandomForest_model.py" -x 2 -d "data/preprocessed.pkl" -i [$i] -n "paper__current__exclude2__ranknorm__dde_only" -f '["gemme_aa_wt_p$"]' ; done

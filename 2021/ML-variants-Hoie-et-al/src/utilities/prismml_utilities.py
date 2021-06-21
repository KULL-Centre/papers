# Load
import sys, os
#sys.path.insert(0, "/lindorffgrp-isilon/maghoi/projects/prismdb__analysis__private/scripts")
import pandas as pd
import numpy as np
import prism_machine_learning as prismml

def predictions_csvs_sanity_check(pm, p1, p0):
    pass_control = True

    #Length checks
    if len(p1) != 1:
        print("p1 length not 1:", len(p1))
        pass_control = False
    if len(pm) != len(p0) + 1:
        print("pm length not p1 + 1:", len(pm), len(p1), len(p0))
        pass_control = False

    #Content check
    if (pm.columns != p1.columns).all() or (pm.columns != p0.columns).all():
        print("pm, p1 or p0 non-matching columns:", pm.columns, p1.columns, p0.columns)
        pass_control = False
    if (pm.iloc[-1] != p1.iloc[-1]).all():
        print("pm, p1 first item not matching:", pm.iloc[-1], p1.iloc[-1])
        pass_control = False
    return(pass_control)


from datetime import datetime, date
from joblib import dump, load
def csv_save_record_performance(model, train_X, train_y, valid_X, valid_y, dfs_proc_mave, dfs_names, valid_ix_list, train_ix_list, stats_df, chosen_features, start_time, RUN_NAME = "RUN_NAME", save_model = False, verbose = 0):

    # Save model with name based on valid proteins
    #outname = str(valid_proteins[0]) + "".join(["_"+ str(n) for n in valid_ix_list]) + ".joblib"
    #outpath = "data/ML_trials/testmodel__14may_2020/" + outname
    #dump(model, outpath)

    # Master CSV vars
    valid_proteins = list(prismml.get_proteins_by_idxs(valid_ix_list, dfs_names).index)
    Pr, Sr, MAE, r2 = prismml.test_performance_continuous(model, valid_X, valid_y, include = "all")
    train_ix_list
    valid_ix_list
    y_pred = model.predict(valid_X)
    y_true = valid_y

    features_uniq = list(pd.unique(valid_X.columns.map(lambda x: x[0:3])))

    all_protein_stats = prismml.test_proteins(model, valid_ix_list, dfs_proc_mave,
                                              dfs_names, stats_df,
                                             include = "all")
    all_protein_Pr, all_protein_Sr = all_protein_stats["Ppred"], all_protein_stats["Spred"]


    # Input vars
    run_name = RUN_NAME + "_"
    DATASET_NAME = RUN_NAME + ".pkl.gz"

    # Base path
    base_path = "data/runs/"
    os.makedirs(base_path, exist_ok=True)

    # Create folder for specific run
    today_string = str(date.today()).replace("-", "_")
    model_outdir = base_path + run_name + today_string + "/"
    os.makedirs(model_outdir, exist_ok=True)
    dataset_path = base_path + DATASET_NAME

    # Save model with name based on valid proteins
    model_name = "trained_model_" + str(valid_proteins[0]) + "".join(["_"+ str(n) for n in valid_ix_list]) + ".joblib"
    model_path = model_outdir + model_name
    if save_model: dump(model, model_path)

    # Load CSV file
    predictions_csv = pd.DataFrame({"run_name": np.nan,
                              "date" : np.nan, "duration": np.nan,
                              "valid_proteins":np.nan, "valid_Pearson": np.nan,
                              "valid_Spearman":np.nan,"valid_MAE": np.nan, "valid_r2":np.nan,
                              "all_Pearson":np.nan, "all_Spearman":np.nan,
                              "chosen_features":np.nan,
                              "train_idxs": np.nan, "valid_idxs": np.nan,
                              "dataset_path": np.nan,
                            "model_path": np.nan,
                            "model_params": np.nan,
                              "y_pred": np.nan, "y_true":np.nan}, index = [0])
    predictions_csv_path = model_outdir + "predictions.csv"

    if not os.path.exists(predictions_csv_path):
        print("Existing predictions.csv not found. Saving to", predictions_csv_path)
        predictions_csv.to_csv(predictions_csv_path, index = None)
    else:
        print("Loading predictions file", predictions_csv_path)
        predictions_csv = pd.read_csv(predictions_csv_path)

    # Add new predictions to temporary DF

    # Duration
    duration = str(datetime.now() - start_time)

    new_predictions_ix = predictions_csv.index.max() + 1
    new_predictions = pd.DataFrame({"run_name": str(RUN_NAME),
                             "date" : str(datetime.now()), "duration": str(duration),
                              "valid_proteins": [valid_proteins], "valid_Pearson": Pr,
                              "valid_Spearman": Sr,"valid_MAE": MAE, "valid_r2": r2,
                              "all_Pearson": [list(all_protein_Pr)], "all_Spearman": [list(all_protein_Sr)],
                              "chosen_features": [list(chosen_features)],
                              "train_idxs": [list(train_ix_list)], "valid_idxs": [list(valid_ix_list)],
                             "dataset_path": str(dataset_path),
                            "model_path": str(model_path),
                            "model_params": str(model.get_params()),
                              "y_pred": [list(y_pred)], "y_true": [list(y_true)]}, index = [new_predictions_ix])

    # Concatenate to old CSV and save
    predictions_csv_merged = pd.concat([predictions_csv, new_predictions])
    if predictions_csvs_sanity_check(predictions_csv_merged, new_predictions, predictions_csv):
        print("Saving csv to", predictions_csv_path)
        predictions_csv_merged.to_csv(predictions_csv_path, index = None)
    else:
        new_predictions_csv_path = predictions_csv_path + "_fail"
        print("CSV merge quality control failed. Saving to", new_predictions_csv_path)
        predictions_csv_merged.to_csv(new_predictions_csv_path, index = None)

    return(True)

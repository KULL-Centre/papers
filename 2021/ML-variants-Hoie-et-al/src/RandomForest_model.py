# Load
import sys, os, ast
sys.path.insert(0, "src/utilities")

import pandas as pd
import numpy as np

import prismml_utilities as prism_utilities
import prism_machine_learning as prismml
import PrismData as PrismData

# Inputs
from argparse import ArgumentParser
parser = ArgumentParser(description="Select training, validation proteins")
parser.add_argument("-n", dest="RUN_NAME", default = "paper__RF__all__1", help="RUN_NAME")
parser.add_argument("-l", dest="LINEAR_REGRESSION", default = False, help="Run linear regression model instead of RF model")
parser.add_argument("-d", dest="DATASET_PATH", default = "data/preprocessed.pkl", help="Dataset path")
parser.add_argument("-s", dest="SUBSET", default = False, help="TREES")
parser.add_argument("-i", dest="VALID_IDXS", default = False, help="VALID_IDXS")
parser.add_argument("-j", dest="TRAIN_IDXS", default = False, help="TRAIN_IDXS")
parser.add_argument("-f", dest="FEATURES", default = "all", help="FEATURES")
parser.add_argument("-t", dest="TREES", default = 150, help="TREES")
parser.add_argument("-x", dest="EXCLUDE_MISSING", default=2, help="Exclude missing Rosetta values (1), or Rosetta and GEMME (2)")
parser.add_argument("-v", dest="VERBOSE", default=0, help="Print error messages")
args = parser.parse_args()
verbose = int(args.VERBOSE)

# Duration
from datetime import datetime
start_time = datetime.now()

LINEAR_FLAG = args.LINEAR_REGRESSION
if LINEAR_FLAG: print("Linear regression model set")

# Generate features
print("\nChosen features:", args.FEATURES)
preprocessed_path = str(args.DATASET_PATH)
print("Loading from preprocessed path:", preprocessed_path)

# Big load script
dfs_proc, dfs_names, dfs_raw = prismml.generate_load_preprocessed_datasets("data/preprocessed/prism*.txt",
                                            normalize_mave_only = True,
                                            preprocessed_path=preprocessed_path)

# Only include rows with present MAVE score
for df in dfs_raw: df.dropna(subset = ["score"], inplace = True)
for df in dfs_proc: df.dropna(subset = ["score"], inplace = True)
for df_raw, df_proc in zip(dfs_raw, dfs_proc):
    df_proc["mave_p0"] = df_raw["mave_p0"]

# Extra pre-proc
#dfs_raw_names = dfs_names.copy()
#dfs_proc_subset = pd.Series(prismml.filter_numeric_fillnans_dfs(dfs_proc))

from copy import deepcopy
dfs_raw_scores = dfs_proc.apply(deepcopy)
for i in range(len(dfs_raw_scores)):
    dfs_raw_scores[i]["score"] = dfs_raw[i]["score"]
#dfs_raw_scores = pd.Series(prismml.filter_numeric_fillnans_dfs(dfs_raw_scores))

# stats df
stats_df_all = prismml.generate_stats_df(dfs_proc, dfs_names)

# Extract experiments only with Rosetta and Gemme values above 0.3, and within 0.2 of each other
if args.SUBSET:
    print("\nExtracting subset ... Only including dataset if ddG OR ddE correlation > 0.3")
    min_threshold = 0.3
    filtered_names = []
    for i, gem, ros in zip(enumerate(stats_df_all.index), stats_df_all["Gemme"], stats_df_all["Rosetta"]):
        if gem > min_threshold or ros > min_threshold: filtered_names.append(i[1])

    def extract_subset(dfs_list, dfs_names, subset):
        return(pd.Series([dfs_list[i] for i in range(len(dfs_names)) if dfs_names[i] in subset]))

    dfs_proc_subset = extract_subset(dfs_proc, dfs_names, filtered_names)
    dfs_raw_subset = extract_subset(dfs_raw, dfs_names, filtered_names)
    dfs_raw_scores_subset = extract_subset(dfs_raw_scores, dfs_names, filtered_names)
    dfs_names_subset = extract_subset(dfs_names, dfs_names, filtered_names)
    if verbose >= 1: print("\nExtracted subset ({}):\n{}".format(str(len(dfs_names_subset)), dfs_names_subset))

else:
    dfs_proc_subset = dfs_proc
    dfs_raw_subset = dfs_raw
    dfs_raw_scores_subset = dfs_raw_scores
    dfs_names_subset = dfs_names

#print("\nExtracted subset (" + len(dfs_names_subset) +"):", dfs_names_subset)

if verbose > 1:
    for x, name in zip([dfs_proc_subset, dfs_raw_scores_subset, dfs_names_subset, dfs_raw_subset], ["dfs_proc_subset", "dfs_raw_scores_subset",  "dfs_names_subset", "dfs_raw_subset"]):
        print("\n" + name, "type:", type(x))
        try: print("Len:", len(x), "First item shape:", x[0].shape, type(x[0]))
        except: print("Len:", len(x), "First item len:", len(x[0]), type(x[0]))

stats_df_subset = prismml.generate_stats_df(dfs_proc_subset, dfs_names_subset)

# Check if selected specific features
if args.FEATURES == "all": chosen_features_re = "gemme_aa_p0|gemme_aa_wt_p$|gemme_M_p0|ros_aa_p0|ros_aa_wt_p$|ros_M_p0|mave_wt|mave_any"
else: chosen_features_re = "|".join(ast.literal_eval(args.FEATURES))
include = dfs_proc_subset[0].columns.str.contains("^score$|" + chosen_features_re)

chosen_features = dfs_proc[0].columns[include]
dfs_proc_subset = pd.Series([df[chosen_features] for df in dfs_proc_subset])
dfs_raw_scores_subset = pd.Series([df[chosen_features] for df in dfs_raw_scores_subset])
print(dfs_proc_subset[0].columns)

######
# ML
# Load training datasets
######

# Choose train idxs subset from input args
if args.VALID_IDXS: valid_ix_list = ast.literal_eval(args.VALID_IDXS)
else: valid_ix_list = [0]
print("valid_ix_list:", valid_ix_list, type(valid_ix_list))

if args.TRAIN_IDXS: train_ix_list = ast.literal_eval(args.TRAIN_IDXS)
else: train_ix_list = list(range(len(dfs_proc_subset)))
print("train_ix_list (input, before filtering):", train_ix_list, type(train_ix_list))


def get_dataset(dataset_proc, dataset_raw, dataset_names, chosen_features_re, valid_ix_list):
    # Setup
    train_ix_list = list(range(len(dfs_names_subset)))

    # Select features
    include = dataset_proc[0].columns.str.contains("^score$|" + chosen_features_re)
    chosen_features = dataset_proc[0].columns[include]

    dataset_proc = pd.Series([df[chosen_features] for df in dataset_proc])
    dataset_raw = pd.Series([df[chosen_features] for df in dataset_raw])
    if verbose: print(dataset_proc[0].columns)

    # Numeric values only
    dataset_proc = prismml.filter_numeric_fillnans_dfs(dataset_proc)

    # Removal of all proteins by name in train_ix_list
    train_names = prismml.get_proteins_by_idxs(train_ix_list, dataset_names)
    valid_names = prismml.get_proteins_by_idxs(valid_ix_list, dataset_names)
    train_protein_idxs = prismml.get_idxs_by_proteins(train_names, dataset_names)

    print("Removing validation proteins from training proteins")
    exclude = np.array(train_names.isin(valid_names), dtype = bool)
    train_ix_list = list(pd.Series(train_ix_list)[~exclude])

    train_names = prismml.get_proteins_by_idxs(train_ix_list, dataset_names)
    train_protein_idxs = prismml.get_idxs_by_proteins(train_names, dataset_names)

    valid_names = prismml.get_proteins_by_idxs(valid_ix_list, dataset_names)
    if verbose: print("Training set:\n", train_names, "\nValidation set:\n", valid_names)


    # training sets. Remove validation proteins from training set
    train_X, train_y, valid_X, valid_y = prismml.generate_train_valid_set_dfs(dataset_proc, valid_ix_list = valid_ix_list,
                                                                     train_ix_list = train_ix_list)
    _, _, valid_X_raw, valid_y_raw  = prismml.generate_train_valid_set_dfs(dataset_raw,
                                                                              valid_ix_list = valid_ix_list,
                                                                          train_ix_list = train_ix_list)

    if verbose > 1: print("Train/valid x/y", len(train_X), len(train_y), len(valid_X), len(valid_y))

    return(train_X, train_y, valid_X, valid_y_raw)

chosen_features_re = "ros_aa_p0|ros_aa_wt_p$|ros_M_p0|gemme_aa_p0|gemme_aa_wt_p$|gemme_M_p0|mave_wt|mave_any"
train_X, train_y, valid_X, valid_y_raw = get_dataset(dfs_proc_subset, dfs_raw_scores_subset, dfs_names_subset, chosen_features_re, valid_ix_list)

if int(args.EXCLUDE_MISSING) >= 1 and "ros_aa_wt_p" in train_X.columns:
    print("Filtering out missing Rosetta values")
    print("Before:", train_X.shape, valid_X.shape)
    missing_train = np.array(train_X["ros_aa_wt_p"] == -100, dtype = bool)
    train_X, train_y  = train_X[~missing_train], train_y[~missing_train]
    missing_val = np.array(valid_X["ros_aa_wt_p"] == -100, dtype = bool)
    valid_X, valid_y_raw = valid_X[~missing_val],  valid_y_raw[~missing_val]
    print("After:", train_X.shape, valid_X.shape)

if int(args.EXCLUDE_MISSING) >= 2 and "gemme_aa_wt_p" in train_X.columns:
    print("Filtering out missing GEMME values")
    print("Before:", train_X.shape, valid_X.shape)
    missing_train = np.array(train_X["gemme_aa_wt_p"] == -100, dtype = bool)
    train_X, train_y  = train_X[~missing_train], train_y[~missing_train]
    missing_val = np.array(valid_X["gemme_aa_wt_p"] == -100, dtype = bool)
    valid_X, valid_y_raw = valid_X[~missing_val], valid_y_raw[~missing_val]
    print("After:", train_X.shape, valid_X.shape)

# Set model
if LINEAR_FLAG:
    print("Linear regression model set")
    from sklearn.linear_model import LinearRegression
    model = LinearRegression()
else:
    print("Random Forest model set")
    from sklearn.ensemble import RandomForestRegressor
    model = RandomForestRegressor(n_estimators = int(args.TREES), max_features = "sqrt", min_samples_leaf=15)

# Train model
model.fit(train_X, train_y)

# Dump model to file
from joblib import dump, load
dump(model, "trained_RF_model.joblib")

# Evaluate
print("Train performance spearman (norm):", prismml.test_performance_continuous(model, train_X, train_y, include = "spearman"))
print("Valid performance spearman (raw):", prismml.test_performance_continuous(model, valid_X, valid_y_raw, include = "spearman"))
print("Valid performance spearman (raw -> recording...):", prismml.test_performance_continuous(model, valid_X, valid_y_raw, include = "spearman"))

# Record performance
exit_value = prism_utilities.csv_save_record_performance(model, train_X, train_y, valid_X, valid_y_raw, dfs_raw_scores_subset, dfs_names_subset,
                                                         valid_ix_list, train_ix_list, stats_df_subset, chosen_features, start_time,
                                                         RUN_NAME = args.RUN_NAME, verbose = 1)

if exit_value != True:
    print("ERROR: csv_save_record_performance exited with error code")

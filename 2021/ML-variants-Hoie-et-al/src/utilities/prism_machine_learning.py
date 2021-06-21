import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob, sys, os, re
from scipy.stats import pearsonr

import matplotlib.pyplot as plt
import seaborn as sns
plt.style.use('seaborn-whitegrid')

#from scripts import *
import prism_machine_learning as prismml

import PrismData
prismparser = PrismData.PrismParser()


aa_order_alphabetical = pd.Series(["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",
           "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
aa_order_alphabetical_index_series = pd.Series(range(len(aa_order_alphabetical)), index = aa_order_alphabetical)

aa_list = list("ARNDCEQGHILKMFPSTWYV")
aa_dict = dict()
for i, aa in enumerate(aa_list):
    aa_dict[aa] = i

aa_dict_alph = dict()
for i, aa in enumerate(aa_order_alphabetical):
    aa_dict_alph[aa] = i


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVR
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsRegressor

def extract_wt_pos_mut_to_df(df, variant_col = "variant"):
    variants = df[variant_col]

    first = variants.apply(lambda x: x[0])
    mid = variants.apply(lambda x: x[1:-1])
    last = variants.apply(lambda x: x[-1])


    df[["wt", "pos", "mut"]] = pd.DataFrame({"wt":first, "pos":mid, "mut":last})
    df['pos']=df["pos"].astype(int)

    # Check whether to exclude data
    exclude = np.zeros(len(df), dtype=bool)
    for i in range(len(df)):
        try: int(mid[i])
        except: exclude[i] = 1
    if sum(exclude) >= 1: print("Skipping", int(sum(exclude)), "rows (data between first and last character not a numerical position)")
    df = df[~exclude]

    return(df)

def generate_mutation_matrix(df, variant_col = "variant", score_col = "score", count = False):
    df = extract_wt_pos_mut_to_df(df)

    array = np.zeros((20,20))
    array[:] = np.nan
    for i, wt, in enumerate(aa_order_alphabetical):
        for i2, mut in enumerate(aa_order_alphabetical):
            #if wt == mut: continue

            wt_list = df["wt"] == wt
            mut_list = df["mut"] == mut
            both = np.logical_and(wt_list, mut_list)

            median_score = np.nanmean(df[score_col][both])
            if count: median_score = sum(both)
            array[i, i2] = median_score
    return(array)


def generate_mutation_matrix_dfs(dfs, protein_names, col_variant = "variant", col_score = "score", plot = False, count = False):
    mutation_arrs = []

    for i in range(len(dfs)):
        print(protein_names[i])

        if not count: mutation_matrix = generate_mutation_matrix(dfs[i])
        else: mutation_matrix = generate_mutation_matrix(dfs[i], count = True)

        mutation_arrs.append(mutation_matrix)

        if plot: plot_mutation_map(mutation_matrix)

    if not count: master_mutation_arr = np.nanmean(mutation_arrs, axis = 0)
    else: master_mutation_arr = np.nansum(mutation_arrs, axis = 0)

    print("Returning np.nanmean merged mutation array, and all individual mutation arrays")
    return(master_mutation_arr, mutation_arrs)

def dataset_x_y(df, y_col = "score"):
    dataset_y = df.iloc[:, 0]
    dataset_x = df.iloc[:, 1:]
    return(dataset_x, dataset_y)

def train_valid_split(df, valid_percent = 0.25):
    split_point = int(len(df) * (1-valid_percent))
    train_df, valid_df = df[0: split_point], df[split_point: ]
    return(train_df, valid_df)

def generate_train_valid_sets(df, valid_percent = 0.25):
    train_df, valid_df = train_valid_split(df, valid_percent)

    train_x, train_y = dataset_x_y(train_df)
    valid_x, valid_y = dataset_x_y(valid_df)

    return(train_x, train_y, valid_x, valid_y)

from sklearn import metrics
from sklearn.metrics import confusion_matrix
from sklearn.metrics import matthews_corrcoef
import scipy as sc
from sklearn.metrics import precision_score
from scipy.stats import spearmanr


from scipy.stats import pearsonr, spearmanr
def return_correlation(test_x, test_y, correlation = "Pearson"):
    exclude_y = np.isnan(test_y)
    exclude_x = np.isnan(test_x)
    exclude_both = np.logical_or(exclude_x, exclude_y)

    #present_fraction = np.round(sum(exclude_both)/sum(exclude_y), 2)
    #if present_fraction < 0.5: print("> 50 % missing:", present_fraction)

    x = np.array(test_x).flatten()[~exclude_both]
    y = np.array(test_y).flatten()[~exclude_both]

    if len(x) < 50 or len(y) < 50:
        return(-100)
    else:
        if correlation == "Pearson": p = np.round(pearsonr(x, y), 3)[0]
        if correlation == "Spearman": p = np.round(spearmanr(x, y), 3)[0]
        return(p)



# Features
def generate_wt_mut_matrix(df):
    arr = np.array(df)
    wt_aas = np.apply_along_axis(lambda x: str(x[0])[0], 1, arr)
    mut_aas = np.apply_along_axis(lambda x: str(x[0][-1:]), 1, arr)
    scores = arr[:, 1]

    fill_df = pd.DataFrame(index = aa_order_alphabetical, columns = aa_order_alphabetical)
    for wt_aa in aa_order_alphabetical:
        wt_aa_pos = wt_aas == wt_aa

        for mut_aa in aa_order_alphabetical:
            mut_aa_index = mut_aas[wt_aa_pos] == mut_aa
            mut_scores = scores[wt_aa_pos][mut_aa_index]
            mut_scores = np.array(mut_scores, dtype = float)
            if len(mut_scores) == 0: continue

            try: median_mut_at_wt = np.nanmedian(mut_scores)
            except:
                print(mut_scores)
                median_mut_at_wt = np.nanmedian(mut_scores)
                median_mut_at_wt = np.nan
            fill_df.loc[wt_aa, mut_aa] = median_mut_at_wt

    fill_df.astype(float)
    return(fill_df)


from sklearn.preprocessing import MinMaxScaler
def find_most_common_score(df, plot = False):
    df_arr = np.array(df).reshape(-1, 1)
    scaler = MinMaxScaler()
    scaler.fit(df_arr)

    # Scale to 0 -1
    df_scaled = scaler.transform(df_arr).flatten()
    values_sorted = np.sort(df_scaled)

    # Generate bin counts
    counts, bins = np.histogram(df_scaled, bins = 50)

    # Find most common bin, calculate median score between previous and next bin
    half_length = int(len(counts)/2)
    half_counts = counts[half_length:]
    ix = np.argmax(half_counts) + half_length
    start, end = np.sum(counts[0:ix]), np.sum(counts[0:ix+1])
    median_value = np.median(values_sorted[start : end])

    if plot:
        print(len(counts))
        print(start, end, ix)
        sns.distplot(df_scaled, bins = 50)
        plt.axvline(median_value, 0, 3, color = "red", lw = 4, ls = "--", alpha = 0.5)
        title_str = str(df.name) + ": " + str(np.round(median_value, 3))
        plt.title(title_str)
    return(median_value)


from scipy.signal import find_peaks
def norm_rightmost_peak(df, plot = False, h_factor = 0.025):
    # Generate bin counts
    counts, bins = np.histogram(df, bins = 50)

    # Find rightmost peak
    h = int(len(df) * h_factor)
    peaks, heights = find_peaks(counts, height=h)

    # Find rightmost peak bin, calculate median score between previous and next bin
    ix = peaks[-1]
    start, end = np.sum(counts[0:ix]), np.sum(counts[0:ix+1])

    values_sorted = np.sort(df)
    median_value = np.median(values_sorted[start : end])

    if plot:
        print(len(counts))
        print(start, end, ix)
        sns.distplot(df, bins = 50)
        plt.axvline(median_value, 0, 3, color = "red", lw = 4, ls = "--", alpha = 0.5)
        title_str = str(df.name) + ": " + str(np.round(median_value, 3))
        plt.title(title_str)
    return(median_value)


from sklearn.preprocessing import QuantileTransformer
def norm_ranknorm(df, n_quantiles = 4):
    quantile_transformer = preprocessing.QuantileTransformer(n_quantiles = n_quantiles)
    transformed = quantile_transformer.fit_transform(df)
    return(transformed)

from sklearn.preprocessing import MinMaxScaler
def norm_minmax_onefeature(df, feature_range = (0, 1)):
    scaler = MinMaxScaler()
    arr = np.array(df)
    if len(arr.shape) >= 2:
        if arr.shape[1] > 1 or len(arr.shape) > 2: print("More than 1 feature detected?", arr.shape)

    arr = arr.reshape(-1, 1)
    return(scaler.fit_transform(arr).flatten())

####
# Features contrinued
###

from sklearn import preprocessing
def norm_quantnorm(df, n_quantiles = 4):
    quantile_transformer = preprocessing.QuantileTransformer(n_quantiles = n_quantiles)
    transformed = quantile_transformer.fit_transform(df)
    return(transformed)


def norm_quantnorm_dfs(dfs_raw,
                       columns = ["score", "gemme_score"],
                      n_quantiles = 4):
    dfs_proc = []
    if columns == "all":
        columns = dfs_raw[0].select_dtypes(include=[np.number]).columns

    for df in dfs_raw:
        df_proc = df.copy()
        df_proc[columns] = norm_quantnorm(df[columns], n_quantiles = n_quantiles)
        dfs_proc.append(df_proc)
    return(dfs_proc)

def valid_SAV(variant):
    return bool(re.search(r'^[A-Z]\d{1,4}[A-Z]$', variant))

def remove_invalid_variants_dfs(dfs, variant_col = "variant"):
    dfs_proc = []

    for df in dfs:
        include = np.ones(len(df), dtype = bool)

        for i, variant in enumerate(df[variant_col]):
            include[i] = valid_SAV(variant)

        dfs_proc.append(df[include])

    return(dfs_proc)


## Test


def generate_wt_mut_matrix(df):
    arr = np.array(df)
    wt_aas = np.apply_along_axis(lambda x: str(x[0])[0], 1, arr)
    mut_aas = np.apply_along_axis(lambda x: str(x[0][-1:]), 1, arr)
    scores = arr[:, 1]

    fill_df = pd.DataFrame(index = aa_order_alphabetical, columns = aa_order_alphabetical)
    for wt_aa in aa_order_alphabetical:
        wt_aa_pos = wt_aas == wt_aa

        for mut_aa in aa_order_alphabetical:
            mut_aa_index = mut_aas[wt_aa_pos] == mut_aa
            mut_scores = scores[wt_aa_pos][mut_aa_index]
            mut_scores = np.array(mut_scores, dtype = float)
            if len(mut_scores) == 0: continue

            try: median_mut_at_wt = np.nanmedian(mut_scores)
            except:
                print(mut_scores)
                median_mut_at_wt = np.nanmedian(mut_scores)
                median_mut_at_wt = np.nan
            fill_df.loc[wt_aa, mut_aa] = median_mut_at_wt

    fill_df.astype(float)
    return(fill_df)

def extract_wt_pos_mut_to_df(df, variant_col = "variant"):
    if variant_col == "index":
        df["variant"] = df.index
        variant_col = "variant"
    variants = df[variant_col]

    first = variants.apply(lambda x: x[0])
    mid = variants.apply(lambda x: x[1:-1])
    last = variants.apply(lambda x: x[-1])


    df[["wt", "pos", "mut"]] = pd.DataFrame({"wt":first, "pos":mid, "mut":last})
    df['pos'] = df["pos"].astype(int)
    return(df)

def create_feature_map(df, score_col = "gemme_score"):
    df = df[["variant", score_col]]
    df = extract_wt_pos_mut_to_df(df)

    #mave_positions = np.array(df["variant"].apply(lambda x: x[1:-1]), dtype=int)
    #mave_variants = np.array(df["variant"].apply(lambda x: x[-1]))
    #mave_sequence_map_order = [aa_order_alphabetical_index_series[aa] for aa in df["wt"]]
    mave_variants_column_index = [aa_order_alphabetical_index_series[aa] for aa in df["mut"]]

    last_residue_pos = np.sort(df["pos"])[-1]

    # Create empty sequence map
    seq_map = np.zeros(shape = (last_residue_pos, 20))
    seq_map[:] = np.nan

    # Fill with variant scores
    for score, pos, aa in zip(df[score_col], df["pos"], mave_variants_column_index):
        seq_map[pos-1, aa] = score

    #residue_means = np.nanmedian(seq_map, axis = 0)
    #position_means = np.nanmedian(seq_map, axis = 1)

    seq_map[np.isnan(seq_map)] = -100
    #position_means[np.isnan(position_means)] = -100

    return(seq_map)

def extract_features(seq_map, df, nullmave = False):
    df = extract_wt_pos_mut_to_df(df)
    #last_residue_pos = np.max(df["pos"])

    #last_residue_pos = np.array(df["variant"].apply(lambda x: x[1:-1]), dtype=int)[-1]
    mave_positions = list(range(np.min(df["pos"]), np.max(df["pos"])))
    mave_residues = np.array(df["wt"])
    mave_mutants = np.array(df["mut"])

    #mave_residues = np.array(df["variant"].apply(lambda x: x[0]))
    #mave_mutants = np.array(df["variant"].apply(lambda x: x[-1]))
    #print(mave_residues[5:10])

    # Mave residues filler
    firstres, lastres = mave_residues[0], mave_residues[-1]
    mave_residues = np.append(mave_residues, [firstres]*3)
    mave_residues = np.append([lastres]*3, mave_residues)
    #print(mave_residues[5:10])

    # Create filler feature map
    filler_map = np.zeros(shape = (seq_map.shape[0]+6, seq_map.shape[1]))
    filler_map[:] = -100
    filler_map[3:-3] = seq_map

    # No filler needed
    residue_means = np.nanmedian(seq_map, axis = 0)

    # Position means
    position_means = np.nanmedian(seq_map, axis = 1)
    filler_position_means = np.append(position_means, [-100]*3)
    filler_position_means = np.append([-100]*3, filler_position_means)

    #Generate WT->Mut matrix
    wt_mut_matrix = np.array(generate_wt_mut_matrix(df), dtype = float)

    # Load global 20x20 avg MAVE-score mutant matrix
    master_mutant_matrix = np.load("data/mut_matrix_alphabetical.npy")

    # Iterate through positions
    outputs = []
    for pos_0 in range(len(seq_map)):
        pos = pos_0 + 3

        for aa in range(len(aa_order_alphabetical)):
            # Numeric index instead of string
            wt = aa_dict_alph[mave_residues[pos]]
            mut = aa

            if not nullmave:
                variants = filler_map[pos].flatten() # Whole positions
                variant_area_x7 = filler_map[pos-3:pos+4, aa] # That aa
                mean_pos_x7 = filler_position_means[pos-3:pos+4] # Mean position score
                mean_aa = residue_means[aa] # Mean for that aa mutation

                wt_any = np.nanmedian(wt_mut_matrix[wt, :])
                #wt_mut = wt_mut_matrix[wt, mut]
                any_mut = np.nanmedian(wt_mut_matrix[:, mut])

                output_vector = np.hstack([variants, variant_area_x7, mean_pos_x7,
                                           mean_aa, wt_any, any_mut])

            # Find WT residue 3 residues back to 3 forwards
            # Features from master mutation matrix, WT -> Mut, WT -> any, any -> Mut
            if nullmave:
                nullmave_wt_to_mut = master_mutant_matrix[wt, aa ]
                nullmave_any_to_mut = np.nanmean(master_mutant_matrix[:, aa ])

                nullmave_wt_pos_x7 = []
                for wt_x in mave_residues[pos-3 : pos+4]:
                    #print(wt_x)
                    mean = np.nanmean(master_mutant_matrix[aa_dict_alph[wt_x], : ])
                    nullmave_wt_pos_x7.append(mean)
                #nullmave_wt_pos_x7 = [] = [ np.nanmean(master_mutant_matrix[aa_dict_alph[wt_x], : ])
                #                  for wt_x in mave_residues[pos-3 : pos+4] ]
                #print(pos, wt_x, nullmave_wt_pos_x7)

                output_vector = np.hstack([nullmave_wt_to_mut, nullmave_any_to_mut, nullmave_wt_pos_x7])

            outputs.append(output_vector)

    return(np.vstack(outputs))

def name_feature_map(feature_map, df, nullmave = False):
    df = extract_wt_pos_mut_to_df(df)
    last_residue_pos = np.sort(df["pos"])[-1]

    # Create dict for WT sequence and positions
    pos_only = pd.unique(np.array(df["variant"].apply(lambda x: x[1:-1])))
    wt_and_pos = pd.unique(np.array(df["variant"].apply(lambda x: x[0:-1])))
    pos_wt_dict = {p : p_wt for p, p_wt in zip(pos_only, wt_and_pos)}

    # Create range from 1 to last mave position
    #last_residue_pos = np.array(df["variant"].apply(lambda x: x[1:-1]), dtype=int)[-1]
    mave_positions = list(range(1, last_residue_pos+1))

    # Create variant names (M1A, M1C ... W200Y) from position 1 to last
    variant_names = []
    for pos in mave_positions:
        for aa in aa_order_alphabetical:
            try: wt_pos = pos_wt_dict[str(pos)]
            except: wt_pos = str(pos)
            variant_names.append(str(wt_pos) + str(aa))

    #position_variants = ["aa_p0_" + aa for aa in aa_order_alphabetical]
    feature_names = (["aa_p0_" + aa for aa in aa_order_alphabetical]
                      + ["aa_wt_p-3", "aa_wt_p-2", "aa_wt_p-1",
                      "aa_wt_p", "aa_wt_p1", "aa_wt_p2",
                      "aa_wt_p3",
                      "M_p-3", "M_p-2", "M_p-1",
                      "M_p0", "M_p1", "M_p2",
                      "M_p3",
                      "mean_aa", "WT_any", "any_to_mut"])

    if nullmave: feature_names = ["nullmave_wt_to_mut", "nullmave_any_to_mut",
                      "nullmave_wt_p-3", "nullmave_wt_p-2", "nullmave_wt_p-1",
                      "nullmave_wt_p", "nullmave_wt_p1", "nullmave_wt_p2",
                      "nullmave_wt_p3"]

    #feature_names = (["aa_p-2", "aa_p-1", "aa_p0", "aa_p1", "aa_p2",
    #                  "M_aa",
    #                    "M_p-2", "M_p-1", "M_p0", "M_p1", "M_p2"]
    #                 + ["aa_p-2_" + aa for aa in aa_order_alphabetical]
    #                 + ["aa_p-1_" + aa for aa in aa_order_alphabetical]
    #                 + ["aa_p0_" + aa for aa in aa_order_alphabetical]
    #                 + ["aa_p1_" + aa for aa in aa_order_alphabetical]
    #                 + ["aa_p2_" + aa for aa in aa_order_alphabetical]
    #                 + ["WT_to_any", "WT_to_mut", "any_to_mut"]
    #                )

    named_feature_map = pd.DataFrame(feature_map, index = variant_names, columns = feature_names)
    return(named_feature_map)

def generate_nullmave_dfs(dfs, dfs_names, variant_col = "variant",
                         master_mutant_matrix_path = "data/mut_matrix_alphabetical.npy"):
    dfs_new = []
    master_mutant_matrix = np.load(master_mutant_matrix_path)
    for df_in, name in zip(dfs, dfs_names):
        df_in.index = df_in[variant_col]
        print(name)
        df = extract_wt_pos_mut_to_df(pd.DataFrame(df_in.index))
        df["mave_wt_to_mut"], df["mave_wt_to_any"], df["mave_any_to_mut"],  = np.nan, np.nan, np.nan

        mave_wt_to_muts = []
        mave_wt_to_anys = []
        mave_any_to_muts = []

        for wt, aa in zip(df["wt"], df["mut"]):
            try:
                wt, aa = aa_order_alphabetical_index_series[wt], aa_order_alphabetical_index_series[aa]
                mave_wt_to_muts.append(master_mutant_matrix[wt, aa ])
                mave_wt_to_anys.append(np.nanmean(master_mutant_matrix[wt, : ]))
                mave_any_to_muts.append(np.nanmean(master_mutant_matrix[:, aa ]))
            except:
                pass

        df["mave_wt_to_mut"], df["mave_wt_to_any"], df["mave_any_to_mut"],  = mave_wt_to_muts, mave_wt_to_anys, mave_any_to_muts

        unique_pos_index = df["pos"].drop_duplicates().index
        df_uniq = df.iloc[unique_pos_index]
        values = df_uniq["mave_wt_to_any"]
        first, last = pd.Series([values.iloc[0]]*3), pd.Series([values.iloc[-1]]*3)
        values_filler = first.append(values).append(last)

        # Iterate thorugh
        df["mave_pm3"], df["mave_pm2"], df["mave_pm1"] = np.nan, np.nan, np.nan
        df["mave_p0"], df["mave_p1"] = np.nan, np.nan
        df["mave_p2"], df["mave_p3"] = np.nan, np.nan
        cols = df.columns[df.columns.str.contains("mave_p")]


        for p in range(len(values)):
            idxs = df[df["pos"] == p+1].index
            row = list(values_filler.iloc[p : p+7])
            df.loc[idxs, cols] = row

        df = df[df.columns[df.columns.str.contains("mave")]]
        df_in.index.names = ["index"]
        df.index = df_in.index
        df = pd.concat([df_in, df], axis = 1)
        dfs_new.append(df)
    return(pd.Series(dfs_new))


### test


def generate_feature_map(df, score_col, prefix):
    seq_map = create_feature_map(df, score_col)
    seq_map = extract_features(seq_map, df)
    seq_map = name_feature_map(seq_map, df)

    seq_map.columns = prefix + "_" + seq_map.columns
    seq_map.insert(0, "variant", seq_map.index)
    return(seq_map)

def generate_feature_map_dfs(dfs, dfs_names):
    merge_dfs = []

    for df, name in zip(dfs, dfs_names):
        print(name)
        df_gemme = generate_feature_map(df, score_col = "gemme_score", prefix = "gemme")
        df_rosetta = generate_feature_map(df, score_col = "Rosetta_ddg_score", prefix = "ros")

        merge_df = pd.merge(df, df_rosetta, on = "variant", how = "inner")
        merge_df = pd.merge(merge_df, df_gemme, on = "variant", how = "inner")
        merge_dfs.append(merge_df)

    return(merge_dfs)

def remove_prefix(s, prefix):
    return s[len(prefix):] if s.startswith(prefix) else s

def remove_multivariants_dfs(dfs):
    r = r'^[ACDEFGHIKLMNPQRSTVWY][0-9]{1,4}[ACDEFGHIKLMNPQRSTVWY]$' # Search for single mutant variants, up to position length 9999
    new_dfs = [df[df["variant"].map(lambda x: bool(re.match(r, x)))] for df in dfs] # True if match, false else
    return(new_dfs)

def extract_specific_columns_dfs(prism_dfs, keep_col_re = "variant|^score|gemmeRosetta|ss|rsa",
                             columns = ["variant","score", "gemme_score", "Rosetta_ddg_score", "ss"]):
    output_dfs = []
    for prism_df in prism_dfs:
        keep_cols = prism_df.columns[prism_df.columns.str.contains(keep_col_re)]
        prism_df = prism_df[keep_cols]
        prism_df.columns = columns
        output_dfs.append(prism_df)
    return(pd.Series(output_dfs))

def generate_load_preprocessed_datasets(load_path_re = "data/preprocessed/prism_merged*.txt",
                                        preprocessed_path = "data/preprocessed.pkl", force_create_new = False, prevent_create_new = True, prevent_save = False, filetype = "prismdata", normalize = True, normalize_mave_only = True, ranknorm = False):

    if os.path.exists(preprocessed_path) and not force_create_new:
        dfs_proc, dfs_names, dfs_raw = pd.read_pickle(preprocessed_path)
        dfs_names = pd.Series([remove_prefix(name, "prism_merged_") for name in dfs_names])
        print("\nSuccessfully loaded features from", str(preprocessed_path))
    elif prevent_create_new != True:
        print("Existing dataset not found or force_create_new in", str(preprocessed_path), "... Creating new with same path")
        filenames = glob.glob(load_path_re)
        dfs_names = pd.Series([os.path.basename(file) for file in filenames])
        dfs_names = pd.Series([remove_prefix(name, "prism_merged_") for name in dfs_names])

        if filetype == "prismdata":
            dfs_raw = pd.Series([prismparser.read(file).dataframe for file in filenames])
        else:
            dfs_raw = pd.Series([pd.read_csv(file, sep = "\t") for file in filenames])

        print(dfs_raw)


        # Remove multivariants, non amino-acid characters
        dfs_raw = remove_multivariants_dfs(dfs_raw)

        # Reorder in position chronological order
        #dfs_raw = order_prism_dfs(dfs_raw)

        # Extract only columns of interest
        print("Extracting columns")
        dfs_raw = extract_specific_columns_dfs(dfs_raw, keep_col_re = "variant|^score|gemme|Rosetta|ss|rsa",
                             columns = ["variant","score", "gemme_score", "Rosetta_ddg_score", "ss", "rsa"])
        print(dfs_raw)

        dfs_raw = generate_nullmave_dfs(dfs_raw, dfs_names)

        # Remove invalid, quant normalize, generate 112 length feature matrix and merge
        #dfs_proc = remove_invalid_variants_dfs(dfs_raw)
        if normalize:
            if normalize_mave_only:
                if ranknorm:
                    print("Rank normalizing MAVE scores")
                    dfs_proc = norm_quantnorm_dfs(dfs_raw, columns = ["score"], n_quantiles = len(dfs_raw))
                else:
                    print("4 quartile normalizing MAVE scores")
                    dfs_proc = norm_quantnorm_dfs(dfs_raw, columns = ["score"])
            elif ranknorm:
                print("Rank normalizing MAVE scores")
                dfs_proc = norm_quantnorm_dfs(dfs_raw, columns = ["score"], n_quantiles = len(dfs_raw))
            else:
                print("4-quantile normalizing entire dataset. This should be set to MAVE only (?)")
                dfs_proc = norm_quantnorm_dfs(dfs_raw, columns = "all")
        else:
            dfs_proc = dfs_raw

        print("Generating feature maps")
        dfs_proc = pd.Series(generate_feature_map_dfs(dfs_proc, dfs_names))

        # Add index
        for df in dfs_proc:
            df.index = df["variant"]
            df = df.drop("variant", axis = 1, inplace = True)

        if not prevent_save:
            print("Saving to", preprocessed_path)
            pd.Series([dfs_proc, dfs_names, dfs_raw]).to_pickle(preprocessed_path)

    else: print("Prevent create new dataset set to", prevent_create_new)

    return(dfs_proc, dfs_names, dfs_raw)

def generate_stats_df(dfs_raw, dfs_names, method = "spearman"):
    df_stats = pd.DataFrame({"Ros_cov":np.nan, "Rosetta": np.nan, "Gemme": np.nan, "RSA":np.nan, "mave_matrix" : np.nan}, index = dfs_names)
    for df, name in zip(dfs_raw, dfs_names):
        corr = np.round(df.corr(method = method).loc["score"], 3)

        df_at_mave_positions = df[~np.isnan(df["score"])]
        ros_cov = np.round(1 - np.sum(np.isnan(df_at_mave_positions["Rosetta_ddg_score"])) / len(df_at_mave_positions), 3)
        df_stats.loc[name, ["Ros_cov", "Rosetta", "Gemme", "RSA", "mave_matrix"]] = ros_cov, corr["Rosetta_ddg_score"], corr["gemme_score"], corr["rsa"], corr["mave_p0"]

    # Correct signs
    df_stats["Rosetta"] = df_stats["Rosetta"] * -1

    print("Returning correlation in", method, "Ros score = -1 * score")
    return(df_stats)

def filter_numeric_fillnans_dfs(dfs, check_col = "score", fill_value = -100):
    dfs_proc = []
    for df in dfs:
        missing = np.isnan(df[check_col])
        df = df[~missing]
        df = df.select_dtypes(include=[np.number])
        df = df.fillna(fill_value)
        dfs_proc.append(df)
    return(pd.Series(dfs_proc))


def remove_fill_nans_dfs(dfs, check_col = "score", fill = "all", fill_value = -100):
    dfs_proc = []

    for df in dfs:
        missing = np.isnan(df[check_col])
        df = df[~missing]

        if fill == "all": df = df.fillna(fill_value)
        dfs_proc.append(df)
    return(dfs_proc)

#####
# ML
#####


def generate_train_valid_set_dfs(dfs, valid_ix_list, train_ix_list, y_col = "score"):
    #train_ix_list =  pd.Series(list(range(0, len(dfs)))).drop(valid_ix_list)

    train = pd.concat(list(dfs[train_ix_list]))
    valid = pd.concat(list(dfs[valid_ix_list]))

    train_X, train_y = train.loc[:, train.columns != y_col], train.loc[:, y_col]
    valid_X, valid_y = valid.loc[:, valid.columns != y_col], valid.loc[:, y_col]
    return(train_X, train_y, valid_X, valid_y)


def dataset_x_y_dfs(dfs, y_col = "score"):
    Xs, ys = [], []
    for df in dfs:
        X, y = df.loc[:, df.columns != y_col], df.loc[:, y_col]
        Xs.append(X), ys.append(y)
    return(Xs, ys)

#def get_proteins_by_name(names, dfs_names):
#    for search_name in names:
#        for protein_name in dfs_names:
#            z = re.search("prism_mave_[0-9]{3}_([A-Za-z0-9]{3,8})_", protein_name)
#            protein = z.group(1)

def get_proteins_by_idxs(idxs, dfs_names,
                         re_pattern = "[0-9]{3}_([A-Za-z0-9-]{3,8})_"):
    protein_list = []

    for ix in idxs:
        try:
            z = re.search(re_pattern, dfs_names.iloc[ix])
            protein = z.group(1)
        except:
            print("Unable to extract", ix, dfs_names.iloc[ix])
        protein_list.append(protein)

    protein_list = pd.Series(protein_list, index = dfs_names[idxs])
    return(protein_list)

def get_idxs_by_proteins(names, dfs_names):
    match_bool = np.zeros(len(dfs_names), dtype=bool)
    for search_name in names:
        idxs_bool = dfs_names.str.contains(search_name)
        match_bool = np.logical_or(match_bool, idxs_bool)
        #print(idxs_bool)

    match_list = list(dfs_names[match_bool].index)
    match_list = pd.Series(match_list, index = dfs_names[match_list])

    return(match_list)

from sklearn.metrics import mean_absolute_error
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import r2_score

def calculate_performance_continuous(y_pred, y_true, verbose = 0,
                                include = "pearson"):
    mae = np.round(mean_absolute_error(y_true, y_pred), 2)
    pearson = np.round(pearsonr(y_true, y_pred), 3)[0]
    spearman = np.round(spearmanr(y_true, y_pred), 3)[0]
    r2 = np.round(r2_score(y_true, y_pred), 3)

    if verbose >= 1: print("Pearson (default), Spearman, MAE, R2")
    if include == "pearson": return(pearson)
    if include == "spearman": return(spearman)
    if include == "all": return(pearson, spearman, mae, r2)

def test_performance_continuous(m, test_x, test_y, verbose = 0,
                                include = "spearman"):
    y_true = test_y
    y_score = m.predict(test_x)

    #check for rf
    #try:
    #    y_score = m.predict_proba(test_x)
    #except:
    #    y_score = m.predict(test_x)

    y_pred = y_score
    if verbose == 1: print("y_score", y_score, "y_pred", y_pred)

    mae = np.round(mean_absolute_error(y_true, y_pred), 3)
    pearson = np.round(pearsonr(y_true, y_pred), 3)[0]
    spearman = np.round(spearmanr(y_true, y_pred), 3)[0]
    r2 = np.round(r2_score(y_true, y_pred), 3)

    if include == "pearson": return(pearson)
    if include == "spearman": return(spearman)
    if include == "all": return(pearson, spearman, mae, r2)


def test_proteins(model, valid_ix_list, dfs_proc, dfs_names, stats_df, include = "valid_proteins"):
    stats_df["Ppred"] = 0
    stats_df["Spred"] = 0
    stats_df["Delta"] = 0

    if include == "all":
        names = dfs_names
        Xs, ys = dataset_x_y_dfs(dfs_proc)
    else:
        names = dfs_names[valid_ix_list]
        Xs, ys = dataset_x_y_dfs(dfs_proc[valid_ix_list])

    # Prefer using normalized scores as Pearson correlation gets fucked by Z-rank norm HSP82
    #_, ys = dataset_x_y_dfs(dfs_raw[valid_ix_list])

    for X, y, name in zip(Xs, ys, names):
        Pr, Sr, _, _ = test_performance_continuous(model, X, y, include = "all")
        stats_df.loc[name, "Ppred"] = Pr
        stats_df.loc[name, "Spred"] = Sr
        stats_df.loc[name, "Delta"] = Pr - stats_df.loc[name, "Gemme"]

    if include != "all": stats_df = stats_df.sort_values(by="Ppred")[::-1]
    return(stats_df)

from sklearn import metrics
import scipy as sc

# binary class

def return_binary_class(train_y, threshold = 0.52):
    binary_classes = np.where(train_y > threshold, 1, 0)
    return(binary_classes)

def test_performance_binary(model, test_X, test_y, return_values = "all"):
    y_hat = model.predict(test_X)
    y_true = test_y

    y_hat = np.where(y_hat < 0.5, 1, 0)
    y_true = np.where(y_true < 0.5, 1, 0)

    #Calculate accuracy
    acc = np.round(metrics.accuracy_score(y_true, y_hat), 3)

    #Calculate confusion_matrix, with true and false positives/negatives
    confusion_m =  metrics.confusion_matrix(y_true, y_hat)
    tn, fp, fn, tp = confusion_m.ravel()

    # Next calculate the false true positive rates and precision
    precision = np.round(metrics.precision_score(y_true, y_hat, average='binary'), 3)
    recall = np.round(metrics.recall_score(y_true, y_hat), 3)

    # Calculate the ROC-AUC and Matthews Correlation Coefficient
    fpr, tpr, thresholds = metrics.roc_curve(y_true, y_hat, pos_label=1)
    roc_auc = np.round(metrics.auc(fpr, tpr), 3)
    prec_rec_auc = np.round(metrics.average_precision_score(y_true, y_hat, pos_label = 1), 3)
    mcc = np.round(metrics.matthews_corrcoef(y_true, y_hat), 3)

    # Calculate pearson correlation beween predicted values (y_hat) and true values (y_true)
    pearson = sc.stats.pearsonr(y_true, y_hat)
    pearson_value, pearson_p_value = np.round(pearson[0], 3), pearson[1]

    print("Accuracy:", acc)
    print("ROC AUC:", roc_auc)
    print("MCC:", mcc)
    print("Pearson w/ p-value:", pearson_value, pearson_p_value)
    #print("\nRecall:", tp, "/", tp+fn, "true positives")
    print("Recall", recall)
    print("Precision", precision)

    print("\nConfusion matrix:")
    print("[[tn, fp]]")
    print(confusion_m)
    print("[[fn tp]]")

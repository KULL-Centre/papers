import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

aa_order_alphabetical = pd.Series(["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",
           "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
aa_order_alphabetical_index_series = pd.Series(range(len(aa_order_alphabetical)), index = aa_order_alphabetical)


from sklearn.metrics import mean_absolute_error
#Note: Changing threshold values from 1.95->2.0, 2.95->3.0, 4.47->4.5
def return_array_stats(df, size = 4, x_name = "rosetta_ddg_score", x_thresholds = [-4.01, 2.0, 3.0, 4.5, 16.01],
                            y_name = "gemme_score", y_thresholds = [-0.01, 0.25, 0.50, 0.75, 1.01],
                      y_true = "y_true", y_pred = "y_pred"):
    """
    Divides dataframe into 16 (size x size) sectors according two feature threshold lists,
    and calculates a given statistical value per sector, filled into an array.
    
    Returns four arrays of shape size x size for plotting with sns.heatmap.
    
    Example use:
    arr_proportion, arr_mae, arr_mae_ratio, arr_density = return_array_stats(df)
    sns.heatmap(arr_proportion)
    """
    
    # Defines empty arrays of shape size x size and fills wit
    # Calculate mean_absolute_error per cell between y_true and y_pred
    # for
    # Thresholds defined from model partial dependence plot of Rosetta vs GEMME contribution to predicted fitness
    
    # Initialize empty 4x4 np.arrays 
    size = 4
    plot_array_proportion = np.full((size, size), np.nan)
    plot_array_mae = np.full((size, size), np.nan)
    plot_array_mae_ratio = np.full((size, size), np.nan)
    plot_array_density = np.full((size, size), np.nan)

    MAE_overall = mean_absolute_error(df[y_true], df[y_pred])

    for x in range(size):
        for y in range(size):
            m1 = np.logical_and(df[y_name] >= y_thresholds[y], df[y_name] < y_thresholds[y+1] )
            m2 = np.logical_and(df[x_name] >= x_thresholds[x], df[x_name] < x_thresholds[x+1] )
            mc = np.logical_and(m1, m2)

            v = df[mc]
            perc_low_fitness = np.sum(v[y_true] <= 0.5) / len(v[y_true])
            mae_pred = mean_absolute_error(v[y_true], v[y_pred])
            mae_ratio = (mae_pred - MAE_overall) / MAE_overall
            n = len(v)

            # Fill arrays for later plotting
            plot_array_proportion[y, x] = perc_low_fitness
            plot_array_mae[y, x] = mae_pred
            plot_array_mae_ratio[y, x] = mae_ratio
            plot_array_density[y, x] = n
            
    print("""Returning 4 arrays:
    plot_array_proportion, plot_array_mae, plot_array_mae_ratio, plot_array_density
    """)
    return(plot_array_proportion, plot_array_mae, plot_array_mae_ratio, plot_array_density)

def plot_array_stats(arr, x_name = "Rosetta ddG", y_name = "GEMME ddE", fmt = ".0%", cmap = "coolwarm"):
    # Plot proportion of low fitness variants
    plt.figure(figsize = (10,10))
    fig = sns.heatmap(arr, cmap = cmap, annot = True, annot_kws = {"size" : 27}, fmt = fmt, linewidths = 5)

    # Figure settings
    plt.title(x_name + " vs " + y_name, size = 20)
    fig.set_facecolor('gray')
    plt.ylabel(y_name)
    plt.xlabel(x_name)
    fig.invert_yaxis()

    # Ticks
    plt.xticks(fig.get_xticks(), [' (-)4.0\n->2.0', '   2.0\n->3.0', '   3.0\n->4.5', '   4.5\n->16.0'], rotation = "horizontal", fontsize = 25)
    plt.yticks(fig.get_yticks(), ["0.00-0.25", "0.25-0.50", "0.50-0.75", "0.75-1.00"], rotation = 45, fontsize = 25)

    # Lines
    plt.axhline(y = 1, c = "black", linewidth = 4, linestyle = "--")
    plt.axhline(y = 2, c = "black", linewidth = 4, linestyle = "--")
    plt.axhline(y = 3, c = "black", linewidth = 4, linestyle = "--")

    plt.axvline(x = 1, c = "black", linewidth = 4, linestyle = "--")
    plt.axvline(x = 2, c = "black", linewidth = 4, linestyle = "--")
    plt.axvline(x = 3, c = "black", linewidth = 4, linestyle = "--")
    
    


def generate_seqmap(df_variant, df_score, verbose = 0):
    arr = np.array(df_variant).reshape(-1, 1)
    index_positions = np.apply_along_axis(lambda x: int(str(x[0])[1:-1]), 1, arr)
    mut_positions = np.apply_along_axis(lambda x: str(x[0][-1:]), 1, arr)
    score_positions = df_score

    # Extract start end
    start = int(np.sort(index_positions)[0])
    end = int(np.sort(index_positions)[-1])
    protein_sequence_positions = pd.Series(list(range(start, end+1)))

    #Plot
    seq_map = np.zeros(shape = (protein_sequence_positions.iloc[-1], 20))
    seq_map[:] = np.nan

    for i, pos, mut, score in zip(range(len(arr)), index_positions, mut_positions, score_positions):
        if mut in "*=~X": continue
        seq_map[pos-1, aa_order_alphabetical_index_series[mut]] = score

    # Remove all empty positions before the first MAVE position
    seq_map = seq_map[start-1 : end+1]

    #Transpose
    seq_map = seq_map.transpose()

    position_means = np.nanmedian(seq_map, axis = 0).reshape(1, -1)
    residue_means = np.nanmedian(seq_map, axis = 1).reshape(-1, 1)

    if verbose >= 1:
        print("Returning seq_map, position_means, residue_means")
        print(seq_map.shape, protein_sequence_positions.shape, residue_means.shape)
    return(seq_map, protein_sequence_positions, position_means, residue_means)

def extract_positions_from_seqmaps(dfs, extract_col, variant_col = "variant", score_col = "y_true", rosetta_col = "rosetta_ddg_score", missing_value = np.nan, required_present_mave = 15, required_present_rosetta = 19):

    data_mave, data_new, choice_list = [], [], []

    # Find and select only positions that have > (15) mave variant scores and at least (19) Rosetta variant scores
    for df in dfs:
        # Create combined masked identifying missing values in MAVE score column and in Rosetta column
        seq_map1 = generate_seqmap(df[variant_col], df[score_col])[0]
        if np.isnan(missing_value): variants_present_pos1 = np.sum(~np.isnan(seq_map1), axis = 0)
        else: variants_present_pos1 = np.sum(~(seq_map1 == missing_value), axis = 0)        
        choices1 = variants_present_pos1 >= required_present_mave

        seq_map2 = generate_seqmap(df[variant_col], df[rosetta_col])[0]
        if np.isnan(missing_value): variants_present_pos2 = np.sum(~np.isnan(seq_map2), axis = 0)
        else: variants_present_pos2 = np.sum(~(seq_map2 == missing_value), axis = 0)
        choices2 = variants_present_pos2 >= required_present_rosetta
        
        choices= np.logical_and(choices1, choices2)

        # Create seqmap from extracted column, but filter to only selected positions
        filtered_seq_maps = []
        for col in [extract_col]:
            # Seq map
            seq_map = generate_seqmap(df[variant_col], df[col])[0]
            filtered_seq_map = seq_map[:, choices].transpose()
            filtered_seq_maps.append(filtered_seq_map)

        data_new.append(filtered_seq_maps[0])
        choice_list.append(choices)
        #data_mave.append(filtered_seq_maps[0])
        #data_new.append(filtered_seq_maps[1])

    data = np.vstack(data_new)
    return(data)

def PCA_analysis(df, n_components = 2):
    df = pd.DataFrame(df)
    
    #Define the number of components we want
    pca = PCA(n_components=n_components)

    #Extract principal component 1 and 2 values for each of our samples
    principalComponents = pca.fit_transform(df)
    principalDf = pd.DataFrame(data = principalComponents)

    #Explained variance
    pca_explained_variance = round(sum(pca.explained_variance_ratio_[0:n_components])*100, 2)
    print("The top", n_components, "principal components explain", pca_explained_variance, " % of the variance")
    
    return(np.array(principalDf))


def gaussian_convolutional_heatmap(x, y, c, conv_std = 4, return_values = False, cmap = "coolwarm_r",
                                   vmin = None, vmax = None, center = None,
                                  title = "MAVE mutational landscape:\nUMAP Gaussian convolution map of protein position features\n"):
    """
    Plots a 2D heatmap by converting features x and y to 2D co-ordinates and coloring by feature c.
    Interpolation of points done using astropy.convolution.Gaussian2DKernel
    
    Example of usage:
    gaussian_convolutional_heatmap(df["Rosetta_ddg_score"], df["gemme_score"], c = df["score"], conv_std = 1, cmap = "RdBu")
    """
    
    df_plot = pd.DataFrame({"x":x, "y":y, "score":c})

    from sklearn.preprocessing import MinMaxScaler
    scaler = MinMaxScaler()
    X = np.array(x).reshape(-1, 1)
    X = scaler.fit_transform(X).flatten()
    X = np.array(X * 1000, dtype = int)

    Y = np.array(y).reshape(-1, 1)
    Y = scaler.fit_transform(Y).flatten()
    Y = np.array(Y * 1000, dtype = int)

    df_plot["X"], df_plot["Y"] = X, Y
    
    arr = np.zeros(shape = (1001, 1001))
    arr[:] = np.nan

    for x, y, c in np.array(df_plot[["X", "Y", "score"]]):
        arr[int(y), int(x)] = c
    arr = arr[::-1]
    
    from astropy.convolution import Gaussian2DKernel, convolve
    kernel = Gaussian2DKernel(x_stddev=conv_std)
    # Fix edge bug from convolution setting values to 0
    # First pads edges with 50, fill-value nan, allows 0 bug from convolution
    # then sets padded area to np.nan
    pad = 50
    arr = np.pad(arr, pad, mode = "constant", constant_values = np.nan)

    # Convolve
    heatmap = convolve(arr, kernel)

    # Reset padded area to np.nan
    mask = np.full(arr.shape, True)
    mask[pad:-pad, pad:-pad] = False
    heatmap[mask] = np.nan
    heatmap = heatmap[~mask].reshape(1001, 1001)

    if return_values:
        return(heatmap)
    else:
        if not vmin: vmin = np.nanmin(heatmap)
        if not vmax: vmax = np.nanmax(heatmap)
        if not center: center = np.nanmean([vmin, vmax])
            
        plt.subplots(figsize = (10,10))
        fig = sns.heatmap(heatmap, cmap = cmap, vmin = vmin, vmax = vmax, center = center, rasterized = True)
        plt.title(title)
        fig.set_facecolor('gray')

        # set x and y axes
        ticks = np.arange(0, len(heatmap)+1, int(len(heatmap) / 10))
        labels = np.array(range(ticks.size))/10
        plt.xticks(ticks, labels, rotation = "horizontal")
        plt.yticks(ticks[: -1], labels[::-1][: -1]) # Exclude first 0.0 to avoid overlap with x-axis labels

        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        return(fig)
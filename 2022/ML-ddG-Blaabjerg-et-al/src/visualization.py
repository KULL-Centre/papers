import os

import Bio.PDB.Polypeptide
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import mean_absolute_error

plt.rcParams["figure.dpi"] = 300
plt.rcParams["figure.figsize"] = [8.0, 8.0]
plt.rcParams.update({"font.size": 14})
import glob
import math
import pickle

import mpl_scatter_density
import ptitprince as pt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.offsetbox import AnchoredText
from matplotlib.ticker import FormatStrFormatter

white_viridis = LinearSegmentedColormap.from_list(
    "white_viridis",
    [
        (0, "#ffffff"),
        (1e-20, "#440053"),
        (0.2, "#404388"),
        (0.4, "#2a788e"),
        (0.6, "#21a784"),
        (0.8, "#78d151"),
        (1, "#fde624"),
    ],
    N=256,
)


def learning_curve_cavity(acc_val_list, acc_train_list, loss_train_list):
    """
    Plots the cavity model learning curve
    """
    fig, ax1 = plt.subplots()
    epochs = np.arange(len(acc_val_list)) + 1
    lns1 = ax1.plot(
        epochs, loss_train_list, label="Train: Cross-entropy loss", color="blue"
    )
    ax1.set_ylabel("Cross-entropy loss")
    ax1.set_xticks([1, 3, 5, 7, 9, 11, 13, 15, 17, 19])
    ax1.set_xlabel("Training epochs")
    ax2 = ax1.twinx()
    lns2 = ax2.plot(epochs, acc_train_list, label="Train: Accuracy", color="green")
    lns3 = ax2.plot(epochs, acc_val_list, label="Val: Accuracy", color="orange")
    ax2.set_ylabel("Accuracy")
    lns = lns1 + lns2 + lns3
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs, loc="center right")
    fig.savefig(
        f"{os.path.dirname(os.getcwd())}/output/cavity_models/cavity_learning_curve.pdf",
        bbox_inches="tight",
    )
    plt.close()


def learning_curve_ds(pearson_r_list, val_loss_list, train_loss_list, model_idx):
    """
    Plots the downstream model learning curve
    """
    fig, ax1 = plt.subplots()
    epochs = np.arange(len(pearson_r_list)) + 1
    lns1 = ax1.plot(
        epochs, pearson_r_list, label="Val: Pearson $\\rho$", color="orange"
    )
    ax1.set_ylabel("Pearson $\\rho$")
    ax1.set_xlabel("Training epochs")
    ax1.set_ylim(0.70, 0.80)
    ax1.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax2 = ax1.twinx()
    lns2 = ax2.plot(epochs, val_loss_list, label="Val: MAE$_{F}$ loss", color="green")
    lns3 = ax2.plot(
        epochs, train_loss_list, label="Train: MAE$_{F}$ loss", color="blue"
    )
    ax2.set_ylabel("MAE$_{F}$ loss")
    ax2.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
    lns = lns1 + lns2 + lns3
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs, loc="center right")
    fig.savefig(
        f"{os.path.dirname(os.getcwd())}/output/ds_models/ds_model_{model_idx}/learning_curve.pdf",
        bbox_inches="tight",
    )
    plt.close()


def learning_curve_ds_with_errorbars():
    """
    Plots the downstream model learning curve with error bars
    """

    model_pearson_all = np.zeros((10, 20))
    model_val_loss_all = np.zeros((10, 20))
    model_train_loss_all = np.zeros((10, 20))

    model_files = glob.glob(f"{os.path.dirname(os.getcwd())}/output/ds_models/**/*.pkl")
    for i, model_file in enumerate(model_files):
        with open(f"{model_file}", "rb") as handle:
            e = pickle.load(handle)
            model_pearson_all[i, :] = e[0]
            model_val_loss_all[i, :] = e[1]
            model_train_loss_all[i, :] = e[2]

    fig, ax1 = plt.subplots()
    epochs = np.arange(model_train_loss_all.shape[1]) + 1
    lns1 = ax1.errorbar(
        epochs,
        np.mean(model_pearson_all, axis=0),
        yerr=np.std(model_pearson_all, axis=0),
        label="Val: Pearson $\\rho$",
        color="orange",
    )
    ax1.set_ylabel("Pearson $\\rho$")
    ax1.set_xlabel("Training epochs")
    ax1.set_ylim(0.70, 0.80)
    ax1.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    ax1.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax2 = ax1.twinx()
    lns2 = ax2.errorbar(
        epochs,
        np.mean(model_val_loss_all, axis=0),
        yerr=np.std(model_val_loss_all, axis=0),
        label="Val: MAE$_{F}$ loss",
        color="green",
    )
    lns3 = ax2.errorbar(
        epochs,
        np.mean(model_train_loss_all, axis=0),
        yerr=np.std(model_train_loss_all, axis=0),
        label="Train: MAE$_{F}$ loss",
        color="blue",
    )
    ax2.set_ylabel("MAE$_{F}$ loss")
    ax2.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
    ax2.set_ylim(0.050, 0.110)
    lns = [lns1] + [lns2] + [lns3]
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs, loc="center right")
    fig.savefig(
        f"{os.path.dirname(os.getcwd())}/output/ds_models/ds_learning_curve_with_errorbars.pdf",
        bbox_inches="tight",
    )
    plt.close()


def hist_plot(df, dataset_key):
    """
    Plots the histogram of ddGs from a dataframe
    """

    # Plot all data
    x = df["score_ml"]
    fig, ax = plt.subplots()
    ax.hist(x)
    ax.set_xlabel("ML \u0394\u0394G")
    ax.set_ylabel("Count")
    ax.set_title(f"{dataset_key} - All")
    fig.savefig(
        f"{os.path.dirname(os.getcwd())}/output/{dataset_key}/total_pred_hist.pdf",
        bbox_inches="tight",
    )
    plt.close()

    # Plot individual proteins
    prot_names = df["pdbid"].unique()
    for pdbid in prot_names:
        # Extract protein data
        df_pdbid = df[df["pdbid"] == pdbid]

        # Plot
        x = df_pdbid["score_ml"]
        fig, ax = plt.subplots()
        ax.hist(x)
        ax.set_xlabel("ML \u0394\u0394G")
        ax.set_ylabel("Count")
        ax.set_title(f"{dataset_key} - {pdbid}")
        fig.savefig(
            f"{os.path.dirname(os.getcwd())}/output/{dataset_key}/{pdbid}_pred_hist.pdf",
            bbox_inches="tight",
        )
        plt.close()


def hist_plot_all(df, dataset_key):
    """
    Plots the histogram of ddGs from a dataframe.
    This function is specific to the large cytosolic data set.
    """

    # Plot all data
    x = df["score_ml"]
    print(len(x))
    fig, ax = plt.subplots()
    ax.hist(x, color="blue", bins=100, weights=[1 / 10 ** 3] * len(x))
    ax.set_xlabel("RaSP \u0394\u0394G [kcal/mol]")
    ax.set_ylabel("Variant count x $10^3$")
    fig.savefig(
        f"{os.path.dirname(os.getcwd())}/output/{dataset_key}/total_pred_hist.pdf",
        bbox_inches="tight",
    )
    plt.close()


def make_scatter(x, y, x_vs_y, dataset_key, pdb):
    """
    Plots scatter a scatter plot from a ddG dataframe
    """

    # Compute statistics
    pearson_r = stats.pearsonr(x, y)[0]
    mae = None
    if dataset_key != "VAMP" and (x_vs_y != "exp_vs_ml" or x_vs_y != "exp_vs_ml"):
        mae = mean_absolute_error(x, y)

    # Make plot
    fig = plt.figure()

    # Add data points
    if len(x) > 400:
        ax = fig.add_subplot(1, 1, 1, projection="scatter_density")
        density = ax.scatter_density(x, y, cmap=white_viridis)
        fig.colorbar(density, label="Number of points per pixel")
    else:
        ax = fig.add_subplot(1, 1, 1)
        ax.scatter(x, y, s=2, c="darkblue")

    # Set xy limits, labels and title
    if dataset_key == "Rosetta":
        xylim = ((-1, 7), (-1, 7))
    elif (
        dataset_key == "Protein_G"
        or dataset_key == "ProTherm"
        or "ProTherm_with_Protein_G"
    ):
        xylim = ((-2, 7), (-2, 7))
    elif dataset_key == "VAMP":
        xylim = ((-0.5, 2), (-2, 7))
    ax.set_xlim(xylim[0][0], xylim[0][1])
    ax.set_ylim(xylim[1][0], xylim[1][1])
    if x_vs_y == "exp_vs_ml":
        ax.set_xlabel("Experimental \u0394\u0394G [kcal/mol]")
        if x_vs_y == "exp_vs_ml":
            ax.set_xlabel("VAMP-seq score [a.u.]")
        ax.set_ylabel("RaSP \u0394\u0394G [kcal/mol]")
    elif x_vs_y == "exp_vs_rosetta":
        ax.set_xlabel("Experimental \u0394\u0394G [kcal/mol]")
        if x_vs_y == "exp_vs_ml":
            ax.set_xlabel("VAMP-seq score [a.u.]")
        ax.set_ylabel("Rosetta \u0394\u0394G [kcal/mol]")
    elif x_vs_y == "rosetta_vs_ml":
        ax.set_xlabel("Rosetta \u0394\u0394G [kcal/mol]")
        ax.set_ylabel("RaSP \u0394\u0394G [kcal/mol]")
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))

    # Add textbox
    text = f"Pearson $\\rho$: {pearson_r:.2f}"
    if mae != None:
        text += f"\nMAE: {mae:.2f}"
    anchored_text = AnchoredText(text, loc="upper left")
    ax.add_artist(anchored_text)

    # Save
    if pdb == None:
        plot_name = "total"
    else:
        plot_name = pdb
    fig.savefig(
        f"{os.path.dirname(os.getcwd())}/output/{dataset_key}/{plot_name}_{x_vs_y}_scatter.pdf",
        bbox_inches="tight",
    )
    plt.close()


def make_scatter_combined(df, dataset_key, x_vs_y, prot_names):
    """
    Plots a scatter plot from a ddG dataframe with multiple proteins
    """

    if dataset_key == "ProTherm_with_Protein_G":

        # Plot each protein
        fig = plt.figure()
        colors = ["red", "fuchsia", "blue", "green", "indigo"]

        for i, pdbid in enumerate(prot_names):
            # Select protein data
            df_pdbid = df[df["pdbid"] == pdbid]
            if x_vs_y == "exp_vs_ml":
                x = df_pdbid["score_exp"]
                y = df_pdbid["score_ml"]
            elif x_vs_y == "exp_vs_rosetta":
                x = df_pdbid["score_exp"]
                y = df_pdbid["score_rosetta"]

            y = y[(x >= -1) & (x <= 7)]
            x = x[(x >= -1) & (x <= 7)]

            # Compute statistics
            pearson_r = stats.pearsonr(x, y)[0]

            # Make plot
            ax = fig.add_subplot(1, 1, 1)
            ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))

            # Add data points
            ax = fig.add_subplot(1, 1, 1)
            ax.scatter(
                x,
                y,
                s=2,
                c=colors[i],
                label=f"{prot_names[i]}, Pearson $\\rho$: {pearson_r:.2f}",
            )

            # Set xy limits, labels and title
            xylim = ((-2, 7), (-2, 7))
            ax.set_xlim(xylim[0][0], xylim[0][1])
            ax.set_ylim(xylim[1][0], xylim[1][1])

        if x_vs_y == "exp_vs_ml":
            ax.set_xlabel("Experimental \u0394\u0394G [kcal/mol]")
            ax.set_ylabel("RaSP \u0394\u0394G [kcal/mol]")
        if x_vs_y == "exp_vs_rosetta":
            ax.set_xlabel("Experimental \u0394\u0394G [kcal/mol]")
            ax.set_ylabel("Rosetta \u0394\u0394G [kcal/mol]")

        # Add legends
        plt.legend(loc="upper left")

        # Save
        plot_name = "combined"
        fig.savefig(
            f"{os.path.dirname(os.getcwd())}/output/{dataset_key}/{dataset_key}_{plot_name}_{x_vs_y}_scatter.pdf",
            bbox_inches="tight",
        )
        plt.close()

    if dataset_key == "VAMP":

        # Map PDB ids to protein names
        name_dict = {
            "5BON": "NUDT15 (PDB: 5BON)",
            "2H11": "TPMT (PDB: 2H11)",
            "1D5R": "PTEN (PDB: 1D5R)",
        }

        # Plot each protein
        for i, pdbid in enumerate(prot_names):
            # Select protein data
            df_pdbid = df[df["pdbid"] == pdbid]
            x = df_pdbid["score_exp"]
            y_ml = df_pdbid["score_ml"]
            y_ml = y_ml[(x >= -1) & (x <= 7)]
            y_rosetta = df_pdbid["score_rosetta"]
            y_rosetta = y_rosetta[(x >= -1) & (x <= 7)]
            x = x[(x >= -1) & (x <= 7)]

            # Compute statistics
            pearson_r_ml = stats.pearsonr(x, y_ml)[0]
            pearson_r_rosetta = stats.pearsonr(x, y_rosetta)[0]

            # Make plot
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))

            # Add data points
            ax.scatter(
                x,
                y_ml,
                s=2,
                c="red",
                label=f"RaSP - Pearson $\\rho$: {pearson_r_ml:.2f}",
            )
            ax.scatter(
                x,
                y_rosetta,
                s=2,
                c="blue",
                label=f"Rosetta - Pearson $\\rho$: {pearson_r_rosetta:.2f}",
            )

            # Set xy limits, labels and title
            xylim = ((-0.5, 2), (-1, 7))
            ax.set_xlim(xylim[0][0], xylim[0][1])
            ax.set_ylim(xylim[1][0], xylim[1][1])

            # Set title
            ax.set_title(name_dict[pdbid])

            # if x_vs_y == "exp_vs_ml":
            #    ax.set_xlabel("VAMP-seq score [a.u.]")
            #    ax.set_ylabel("RaSP \u0394\u0394G [kcal/mol]")
            ax.set_xlabel("VAMP-seq score [a.u.]")
            ax.set_ylabel("\u0394\u0394G [kcal/mol]")

            # Add legends
            legend = plt.legend(loc="upper left")

            # Make the font larger for this particular plot
            for item in (
                [ax.title, ax.xaxis.label, ax.yaxis.label]
                + ax.get_xticklabels()
                + ax.get_yticklabels()
                + legend.get_texts()
            ):
                item.set_fontsize(22)

            # Save
            plot_name = "combined"
            fig.savefig(
                f"{os.path.dirname(os.getcwd())}/output/{dataset_key}/{dataset_key}_{plot_name}_{x_vs_y}_scatter_{pdbid}.pdf",
                bbox_inches="tight",
            )
            plt.close()


def scatter_plots(df, dataset_key):
    """
    Function for plotting multiple scatter plots depending on the dataset
    """

    # Create plot: ML vs Rosetta - All
    x = df["score_rosetta"]
    y = df["score_ml"]
    y = y[(x >= -1) & (x <= 7)]
    x = x[(x >= -1) & (x <= 7)]
    make_scatter(x, y, "rosetta_vs_ml", dataset_key, None)

    # Create plot: ML vs Rosetta - Individual
    prot_names = df["pdbid"].unique()
    for pdbid in prot_names:
        # Extract protein data
        df_pdbid = df[df["pdbid"] == pdbid]
        x = df_pdbid["score_rosetta"]
        y = df_pdbid["score_ml"]
        y = y[(x >= -1) & (x <= 7)]
        x = x[(x >= -1) & (x <= 7)]
        make_scatter(x, y, "rosetta_vs_ml", dataset_key, pdbid)

    if dataset_key != "Rosetta" and dataset_key != "Xstal_vs_AF2":
        # Create plot: ML vs Exp - All
        x = df["score_exp"]
        y = df["score_ml"]
        make_scatter(x, y, "exp_vs_ml", dataset_key, None)

        # Create plot: ML vs Exp - Individual
        prot_names = df["pdbid"].unique()
        for pdbid in prot_names:
            # Extract protein data
            df_pdbid = df[df["pdbid"] == pdbid]
            x = df_pdbid["score_exp"]
            y = df_pdbid["score_ml"]
            make_scatter(x, y, "exp_vs_ml", dataset_key, pdbid)

        # Create plot: Rosetta vs Exp - All
        x = df["score_exp"]
        y = df["score_rosetta"]
        make_scatter(x, y, "exp_vs_rosetta", dataset_key, None)

        # Create plot: Rosetta vs Exp - Individual
        prot_names = df["pdbid"].unique()
        for pdbid in prot_names:
            # Extract protein data
            df_pdbid = df[df["pdbid"] == pdbid]
            x = df_pdbid["score_exp"]
            y = df_pdbid["score_rosetta"]
            make_scatter(x, y, "exp_vs_rosetta", dataset_key, pdbid)

    if dataset_key == "ProTherm_with_Protein_G":
        # Create plot: ML vs Exp - All --> Individual colors
        prot_names = df["pdbid"].unique()
        make_scatter_combined(df, dataset_key, "exp_vs_ml", prot_names)
        make_scatter_combined(df, dataset_key, "exp_vs_rosetta", prot_names)

    if dataset_key == "VAMP":
        # Create plot: ML vs VAMP & Rosetta vs VAMP --> Individual colors
        prot_names = df["pdbid"].unique()
        make_scatter_combined(df, dataset_key, "exp_vs_ml", prot_names)


def loss_analysis(df, dataset_key):
    """
    Plots three 20x20 matrices:
       1) Loss per amino acid
       2) ddG variance per amino acid
       3) Residue count in high-loss region per amino acid
    """

    # Order aa by chemistry (aromatic, hydrophobic, hydroxylic, polar, negatively charged)
    aa_list = [
        "W",
        "Y",
        "F",
        "L",
        "M",
        "I",
        "V",
        "G",
        "A",
        "S",
        "T",
        "C",
        "N",
        "H",
        "Q",
        "K",
        "R",
        "D",
        "E",
        "P",
    ]

    # Extract wt and mt
    df["wt"] = df["variant"].str.extract("(^[a-zA-Z])")
    df["mt"] = df["variant"].str.extract("([a-zA-Z]$)")

    # Compute MAEs
    df["mae_ml_rosetta"] = (df["score_ml"] - df["score_rosetta"]).abs()
    if dataset_key != "Rosetta" and dataset_key != "Xstal_vs_AF2":
        df["mae_ml_exp"] = (df["score_ml"] - df["score_exp"]).abs()

    # Divide into high- and low-loss
    df_highloss_rosetta = df[df["mae_ml_rosetta"] > 2]
    df_lowloss_rosetta = df[df["mae_ml_rosetta"] < 2]
    if dataset_key != "Rosetta":
        df_highloss_exp = df[df["mae_ml_rosetta"] > 2]
        df_highloss_exp = df[df["mae_ml_rosetta"] < 2]

    # Compute statistics per mutation
    mat_loss_ml_rosetta = np.zeros((20, 20))
    mat_var_rosetta = np.zeros((20, 20))
    mat_highloss_count_rosetta = np.zeros((20, 20))
    if dataset_key != "Rosetta":
        mat_loss_ml_exp = np.zeros((20, 20))
        mat_var_exp = np.zeros((20, 20))
        mat_highloss_count_exp = np.zeros((20, 20))

    for i, res_wt in enumerate(aa_list):
        for j, res_mt in enumerate(aa_list):
            mat_loss_ml_rosetta[i, j] = df["mae_ml_rosetta"][
                (df["wt"] == res_wt) & (df["mt"] == res_mt)
            ].mean()
            mat_var_rosetta[i, j] = df["score_rosetta_fermi"][
                (df["wt"] == res_wt) & (df["mt"] == res_mt)
            ].var()
            if len(df[(df["wt"] == res_wt) & (df["mt"] == res_mt)]) > 0:
                mat_highloss_count_rosetta[i, j] = len(
                    df_highloss_rosetta[
                        (df_highloss_rosetta["wt"] == res_wt)
                        & (df_highloss_rosetta["mt"] == res_mt)
                    ]
                ) / len(df[(df["wt"] == res_wt) & (df["mt"] == res_mt)])
            if dataset_key != "Rosetta" and dataset_key != "Xstal_vs_AF2":
                mat_loss_ml_exp[i, j] = df["mae_ml_exp"][
                    (df["wt"] == res_wt) & (df["mt"] == res_mt)
                ].mean()
                mat_var_exp[i, j] = df["score_exp_fermi"][
                    (df["wt"] == res_wt) & (df["mt"] == res_mt)
                ].var()
                if len(df[(df["wt"] == res_wt) & (df["mt"] == res_mt)]) > 0:
                    mat_highloss_count_exp[i, j] = len(
                        df_highloss_exp[
                            (df_highloss_exp["wt"] == res_wt)
                            & (df_highloss_exp["mt"] == res_mt)
                        ]
                    ) / len(df[(df["wt"] == res_wt) & (df["mt"] == res_mt)])

    # Make loss, var and count plots
    make_heatmap(dataset_key, mat_loss_ml_rosetta, aa_list, "loss_ml_rosetta")
    make_heatmap(dataset_key, mat_var_rosetta, aa_list, "var_rosetta")
    make_heatmap(
        dataset_key, mat_highloss_count_rosetta, aa_list, "highloss_count_rosetta"
    )
    if dataset_key != "Rosetta" and dataset_key != "Xstal_vs_AF2":
        make_heatmap(dataset_key, mat_loss_ml_exp, aa_list, "loss_ml_exp")
        make_heatmap(dataset_key, mat_var_exp, aa_list, "var_exp")
        make_heatmap(dataset_key, mat_highloss_count_exp, aa_list, "highloss_count_exp")


def make_heatmap(dataset_key, mat, aa_list, description):
    fig, ax = plt.subplots()
    fig.set_size_inches(10.0, 8.0)
    df_heat = pd.DataFrame(mat, columns=aa_list, index=aa_list)
    heatmap_loss = plt.pcolor(df_heat)
    cbar = plt.colorbar(heatmap_loss)
    cbar.set_label("Mean absolute error [kcal/mol]")
    ax.set_yticks(np.arange(0.5, len(df_heat.index), 1))
    ax.set_yticklabels(aa_list)
    ax.set_xticks(np.arange(0.5, len(df_heat.columns), 1))
    ax.set_xticklabels(aa_list)
    ax.set_xlabel("Variant amino acid")
    ax.set_ylabel("Wild-type amino acid")
    fig.savefig(
        f"{os.path.dirname(os.getcwd())}/output/{dataset_key}/heatmap_{description}.pdf",
        bbox_inches="tight",
    )
    plt.close()


def homology_plot(df, dataset_key):
    """
    Plots correlation coefficients to experimental ddGs based on homologous structures.
    Lines for both Rosetta and RaSP.
    """

    # Make new columns
    df["pdb"] = df["pdbid"].str[:4]
    df["bin"] = df["pdbid"].str[8:9].astype("int64")

    pdb_list = df["pdb"].unique()
    bin_list = df["bin"].unique()

    # Plot
    fig, ax_list = plt.subplots(nrows=len(pdb_list), ncols=1, sharex=True)
    fig.subplots_adjust(wspace=0.2)
    fig.tight_layout(pad=2, w_pad=5, h_pad=1)
    fig.text(0.5, -0.001, "Sequence identity", ha="center")
    fig.text(-0.0005, 0.5, "Pearson $\\rho$", va="center", rotation="vertical")

    # Define protein names
    name_dict = {
        "1LZ1": "Lysozyme",
        "2CI2": "Chymotrypsin inhibitor",
        "2RN2": "RNAse",
        "1BVC": "Myoglobin",
    }

    # Define vertical line and seq id positions
    vline_list = [100, 85, 70, 55, 40, 25]

    x_dict = {
        "1LZ1": [100, 98, 78, 62, 47, 33],
        "2CI2": [100, 91, 82, 40, 34],
        "2RN2": [100, 90, 72, 67, 43, 30],
        "1BVC": [100, 92, 81, 65, 43, 31],
    }

    # Lopp over proteins
    for i, pdb in enumerate(pdb_list):

        df_pdb = df[df["pdb"] == pdb]
        df_pdb = df_pdb.sort_values("bin")
        bin_list = df_pdb["bin"].unique()
        ax = ax_list[i]

        # Define x positions
        x = x_dict[pdb]

        # Loop over data points
        corr_ml = []
        corr_rosetta = []
        for i, binid in enumerate(bin_list):
            df_bin = df_pdb[df_pdb["bin"] == binid]
            df_bin_mean = df_bin.groupby("variant").mean()
            if len(df_bin_mean) != 0:
                corr_ml.append(
                    stats.pearsonr(df_bin_mean["score_ml"], df_bin_mean["score_exp"])[0]
                )
                corr_rosetta.append(
                    stats.pearsonr(
                        df_bin_mean["score_rosetta"], df_bin_mean["score_exp"]
                    )[0]
                )

        # Make scatter plot
        l1 = ax.plot(x, corr_rosetta, label="Rosetta", color="blue", marker=".")
        l2 = ax.plot(x, corr_ml, label="RaSP", color="orange", marker=".")
        anchored_text = AnchoredText(
            name_dict[pdb], loc="upper right", prop={"size": 8}
        )
        ax.add_artist(anchored_text)
        ax.set_xticks([92.5, 77.5, 62.5, 47.5, 32.5])
        ax.set_yticks([1.00, 0.75, 0.5, 0.25])
        ax.set_yticklabels(["1.00", "0.75", "0.50", "0.25"], fontsize="8")
        lns = l1 + l2
        labels = [l.get_label() for l in lns]

        # Define dashed vertical lines
        for vline in vline_list:
            ax.axvline(vline, color="lightgrey", linestyle="--")

    # Define x axis text and save fig
    ax.invert_xaxis()
    ax.set_xticklabels(["bin1", "bin2", "bin3", "bin4", "bin5"], fontsize="9")
    ax.text(vline_list[0] + 2.5, 0.05, "100%", color="dimgrey", fontsize="12")
    ax.text(vline_list[1] + 2.5, 0.05, "85%", color="dimgrey", fontsize="12")
    ax.text(vline_list[2] + 2.5, 0.05, "70%", color="dimgrey", fontsize="12")
    ax.text(vline_list[3] + 2.5, 0.05, "55%", color="dimgrey", fontsize="12")
    ax.text(vline_list[4] + 2.5, 0.05, "40%", color="dimgrey", fontsize="12")
    ax.text(vline_list[5] + 2.5, 0.05, "25%", color="dimgrey", fontsize="12")
    fig.legend(lns, labels, loc="upper center", ncol=2, prop={"size": 12})
    fig.savefig(
        f"{os.path.dirname(os.getcwd())}/output/{dataset_key}/homology_plot.pdf",
        bbox_inches="tight",
    )


def gnomAD_category(df):
    if df > 10 ** -2:
        return "gnomAD AF > 1e-2"

    elif df < 10 ** -2 and df > 10 ** -4:
        return "1e-2 > gnomAD AF > 1e-4"

    elif df < 10 ** -4:
        return "gnomAD AF < 1e-4"
    else:
        return np.nan


def category_back(df):
    if df != np.nan:
        return "Back"
    else:
        return np.nan


def intvl_95(np_array):
    d = []
    N = len(np_array)
    p = math.floor(0.95 * N)

    for i in range(N - p):
        d_i = np_array[i + p] - np_array[i]
        d.append(d_i)

    i = np.argmin(d)
    j = i + p

    r = (np_array[i], np_array[j])
    return np.around(r, decimals=2)


def plot_gnomad_clinvar():
    # Load data
    df = pd.read_csv(
        f"{os.path.dirname(os.getcwd())}/data/test/GnomAD_ClinVar/df_rasp_gnomad_clinvar.csv"
    )
    df["gnomAD_category"] = df.apply(lambda x: gnomAD_category(x["gnomad;AF"]), axis=1)
    df["category_back"] = df.apply(lambda x: category_back(x["RaSP;score_ml"]), axis=1)

    # Initialize
    y_RaSP = "RaSP;score_ml"
    x_back = "category_back"
    x_clinvar = "clinvar;Clinvar_signifiance"
    x_gnomad = "gnomAD_category"
    xlim = [-1, 7]

    # Data arrays
    benign_arr = np.sort(
        np.array(
            df["RaSP;score_ml"].values[
                np.where(df["clinvar;Clinvar_signifiance"].values == "benign")
            ]
        )
    )
    pathogenic_arr = np.sort(
        np.array(
            df["RaSP;score_ml"].values[
                np.where(df["clinvar;Clinvar_signifiance"].values == "pathogenic")
            ]
        )
    )

    gnomAD0_arr = np.sort(
        np.array(
            df["RaSP;score_ml"].values[
                np.where(df["gnomAD_category"].values == "gnomAD AF > 1e-2")
            ]
        )
    )
    gnomAD1_arr = np.sort(
        np.array(
            df["RaSP;score_ml"].values[
                np.where(df["gnomAD_category"].values == "1e-2 > gnomAD AF > 1e-4")
            ]
        )
    )
    gnomAD2_arr = np.sort(
        np.array(
            df["RaSP;score_ml"].values[
                np.where(df["gnomAD_category"].values == "gnomAD AF < 1e-4")
            ]
        )
    )

    # Plotting Raincloud plots of different distributions
    # Setting plot style
    sns.set(style="whitegrid", font_scale=1.5)

    # Setting up subplot figure
    fig, axs = plt.subplots(nrows=6, ncols=1, figsize=(10, 20))

    # Background RaSP distribution plotted as a raincloud
    pt.RainCloud(
        x=x_back,
        y=y_RaSP,
        data=df,
        order=["Back"],
        orient="h",
        palette=sns.color_palette("colorblind")[7:],
        bw=0.1,
        alpha=0.65,
        ax=axs[0],
    )

    # ClinVar 'benign' classified RaSP distribution plotted as a raincloud
    pt.RainCloud(
        x=x_clinvar,
        y=y_RaSP,
        data=df,
        order=["benign"],
        orient="h",
        palette=sns.color_palette("colorblind")[0:],
        bw=0.1,
        alpha=0.65,
        ax=axs[1],
    )

    cat_benign_median = np.nanmedian(
        df["RaSP;score_ml"].values[
            np.where(df["clinvar;Clinvar_signifiance"].values == "benign")
        ]
    )
    i_benign, j_benign = intvl_95(benign_arr)

    axs[1].annotate(
        "Median: "
        + str(np.around(cat_benign_median, decimals=2))
        + " "
        + "["
        + str(i_benign)
        + ";"
        + str(j_benign)
        + "]",
        (-0.5, 0.25),
        fontsize=11,
        color=sns.color_palette("colorblind")[0],
    )

    # ClinVar 'pathogenic' classified RaSP distribution plotted as a raincloud
    pt.RainCloud(
        x=x_clinvar,
        y=y_RaSP,
        data=df,
        order=["pathogenic"],
        orient="h",
        palette=sns.color_palette("colorblind")[1:],
        bw=0.1,
        alpha=0.65,
        ax=axs[2],
    )

    cat_patho_median = np.nanmedian(
        df["RaSP;score_ml"].values[
            np.where(df["clinvar;Clinvar_signifiance"].values == "pathogenic")
        ]
    )
    i_patho, j_patho = intvl_95(pathogenic_arr)

    axs[2].annotate(
        "Median: "
        + str(np.around(cat_patho_median, decimals=2))
        + " "
        + "["
        + str(i_patho)
        + ";"
        + str(j_patho)
        + "]",
        (-0.5, 0.25),
        fontsize=11,
        color=sns.color_palette("colorblind")[1],
    )

    # RaSP distribution for gnomAD AF_tot > 1e-2 plotted as a raincloud
    pt.RainCloud(
        x=x_gnomad,
        y=y_RaSP,
        data=df,
        order=["gnomAD AF > 1e-2"],
        palette=sns.color_palette("colorblind")[2:],
        orient="h",
        bw=0.1,
        alpha=0.65,
        ax=axs[3],
    )

    cat_0_median = np.nanmedian(
        df["RaSP;score_ml"].values[
            np.where(df["gnomAD_category"].values == "gnomAD AF > 1e-2")
        ]
    )
    i_cat0, j_cat0 = intvl_95(gnomAD0_arr)

    axs[3].annotate(
        "Median: "
        + str(np.around(cat_0_median, decimals=2))
        + " "
        + "["
        + str(i_cat0)
        + ";"
        + str(j_cat0)
        + "]",
        (-0.5, 0.25),
        fontsize=11,
        color=sns.color_palette("colorblind")[2],
    )

    # RaSP distribution for 1e-2 > gnomAD AF_tot > 1e-4 plotted as a raincloud
    pt.RainCloud(
        x=x_gnomad,
        y=y_RaSP,
        data=df,
        order=["1e-2 > gnomAD AF > 1e-4"],
        palette=sns.color_palette("colorblind")[3:],
        orient="h",
        bw=0.1,
        alpha=0.65,
        ax=axs[4],
    )

    cat_1_median = np.nanmedian(
        df["RaSP;score_ml"].values[
            np.where(df["gnomAD_category"].values == "1e-2 > gnomAD AF > 1e-4")
        ]
    )
    i_cat1, j_cat1 = intvl_95(gnomAD1_arr)

    axs[4].annotate(
        "Median: "
        + str(np.around(cat_1_median, decimals=2))
        + " "
        + "["
        + str(i_cat1)
        + ";"
        + str(j_cat1)
        + "]",
        (-0.5, 0.25),
        fontsize=11,
        color=sns.color_palette("colorblind")[3],
    )

    # RaSP distribution for gnomAD AF_tot < 1e-4 plotted as a raincloud
    pt.RainCloud(
        x=x_gnomad,
        y=y_RaSP,
        data=df,
        order=["gnomAD AF < 1e-4"],
        palette=sns.color_palette("colorblind")[4:],
        orient="h",
        bw=0.1,
        alpha=0.65,
        ax=axs[5],
    )

    cat_2_median = np.nanmedian(
        df["RaSP;score_ml"].values[
            np.where(df["gnomAD_category"].values == "gnomAD AF < 1e-4")
        ]
    )
    i_cat2, j_cat2 = intvl_95(gnomAD2_arr)

    axs[5].annotate(
        "Median: "
        + str(np.around(cat_2_median, decimals=2))
        + " "
        + "["
        + str(i_cat2)
        + ";"
        + str(j_cat2)
        + "]",
        (-0.5, 0.25),
        fontsize=11,
        color=sns.color_palette("colorblind")[4],
    )

    # Global setting for all subplots
    for ax in axs:
        ax.set_xlim(xlim)
        ax.set_xlabel("RaSP \u0394\u0394G [kcal/mol]", fontsize=12, labelpad=15)
        ax.set(ylabel=None)
        ax.tick_params(labelsize=12)

    fig.align_labels()
    fig.savefig(
        f"{os.path.dirname(os.getcwd())}/output/GnomAD_ClinVar/rasp_gnomad_clinvar.pdf",
        bbox_inches="tight",
    )

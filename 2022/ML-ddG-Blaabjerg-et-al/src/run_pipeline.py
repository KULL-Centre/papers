import glob
import os
import pathlib
import random
import subprocess
import sys

import numpy as np
import pandas as pd
import torch
from Bio.PDB.Polypeptide import index_to_one
from torch.utils.data import DataLoader, Dataset

from cavity_model import (
    CavityModel,
    DownstreamModel,
    ResidueEnvironment,
    ResidueEnvironmentsDataset,
)
from helpers import (
    compute_pdb_combo_corrs,
    ds_pred,
    ds_train_val,
    fermi_transform,
    get_ddg_dataloader,
    init_lin_weights,
    inverse_fermi_transform,
    populate_dfs_with_resenvs,
    remove_disulfides,
    train_loop,
    train_val_split_cavity,
    train_val_split_ds,
)
from visualization import (
    hist_plot_all,
    homology_plot,
    learning_curve_ds_with_errorbars,
    loss_analysis,
    plot_gnomad_clinvar,
    scatter_plots,
)

# Set fixed seed
seed = 0
torch.manual_seed(seed)
np.random.seed(seed)
random.seed(seed)
torch.cuda.manual_seed_all(seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
torch.backends.cudnn.enabled = False

# Cavity constants
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
TRAIN_VAL_SPLIT_CAVITY = 0.9
BATCH_SIZE_CAVITY = 100
SHUFFLE_PDBS_CAVITY = True
LEARNING_RATE_CAVITY = 3e-4
EPOCHS_CAVITY = 20
PATIENCE_CUTOFF = 6

# Downstream constants
BATCH_SIZE_DS = 40
LEARNING_RATE_DS = 5e-4
EPOCHS_DS = 20
NUM_ENSEMBLE = 10


def main():
    # Pre-process all protein structures
    print(f"Pre-processing PDBs ...")
    pdb_dirs = [
        f"{os.path.dirname(os.getcwd())}/data/train/cavity/structure/",
        f"{os.path.dirname(os.getcwd())}/data/train/downstream/structure/",
        f"{os.path.dirname(os.getcwd())}/data/test/ProTherm/homology_models/structure/",
        f"{os.path.dirname(os.getcwd())}/data/test/ProTherm/targets/structure/",
        f"{os.path.dirname(os.getcwd())}/data/test/Protein_G/structure/",
        f"{os.path.dirname(os.getcwd())}/data/test/Rosetta_10/structure/",
        f"{os.path.dirname(os.getcwd())}/data/test/VAMP/structure/",
        f"{os.path.dirname(os.getcwd())}/data/test/Xstal_vs_AF2/structure/",
    ]
    for pdb_dir in pdb_dirs:
        subprocess.call(
            [
                f"{os.path.dirname(os.getcwd())}/src/pdb_parser_scripts/parse_pdbs_pred.sh",
                str(pdb_dir),
            ]
        )
    print("Pre-processing finished.")

    # Load parsed PISCES PDBs and perform train/val split
    pdb_filenames_cavity = sorted(
        glob.glob(
            f"{os.path.dirname(os.getcwd())}/data/train/cavity/structure/parsed/*coord*"
        )
    )

    if SHUFFLE_PDBS_CAVITY:
        random.shuffle(pdb_filenames_cavity)

    (
        dataloader_train_cavity,
        dataset_train_cavity,
        dataloader_val_cavity,
        dataset_val_cavity,
    ) = train_val_split_cavity(
        pdb_filenames_cavity, TRAIN_VAL_SPLIT_CAVITY, BATCH_SIZE_CAVITY, DEVICE
    )

    # Define cavity model
    cavity_model_net = CavityModel(get_latent=False).to(DEVICE)
    loss_cavity = torch.nn.CrossEntropyLoss()
    optimizer_cavity = torch.optim.Adam(
        cavity_model_net.parameters(), lr=LEARNING_RATE_CAVITY
    )

    # Train cavity model
    print("Starting cavity model training")
    best_cavity_model_path = train_loop(
        dataloader_train_cavity,
        dataloader_val_cavity,
        cavity_model_net,
        loss_cavity,
        optimizer_cavity,
        EPOCHS_CAVITY,
        PATIENCE_CUTOFF,
    )
    print("Finished cavity model training")

    # Create temporary residue environment datasets to more easily match ddG data
    pdb_filenames_ds = sorted(
        glob.glob(
            f"{os.path.dirname(os.getcwd())}/data/train/downstream/structure/parsed/*coord*"
        )
    )
    dataset = ResidueEnvironmentsDataset(pdb_filenames_ds, transformer=None)
    resenv_dataset = {}
    for resenv in dataset:
        key = (
            f"{resenv.pdb_id}{resenv.chain_id}_{resenv.pdb_residue_number}"
            f"{index_to_one(resenv.restype_index)}"
        )
        resenv_dataset[key] = resenv

    # Load ddG data to dataframe
    df_ddg = pd.read_csv(
        f"{os.path.dirname(os.getcwd())}/data/train/downstream/ddG_Rosetta/ddg.csv"
    )

    # Populate dataframes with wt ResidueEnvironment objects and wt and mt restype indices
    populate_dfs_with_resenvs(df_ddg, resenv_dataset)

    # Remove disulfides
    n_df_start = len(df_ddg)
    df_ddg = remove_disulfides(df_ddg)
    print(
        f"{n_df_start-len(df_ddg)} data points removed due to disulfides in training and validation data."
    )

    # Do Fermi transform
    df_ddg["score_fermi"] = df_ddg["score"].apply(fermi_transform)

    # Checkpoint - save
    df_ddg.to_pickle(f"{os.path.dirname(os.getcwd())}/output/ddg_train_val.pkl")
    f = open(
        f"{os.path.dirname(os.getcwd())}/output/cavity_models/best_model_path.txt", "w"
    )
    f.write(best_cavity_model_path)
    f.close()

    # Checkpoint - load
    df_ddg = pd.read_pickle(f"{os.path.dirname(os.getcwd())}/output/ddg_train_val.pkl")
    best_cavity_model_path = (
        open(
            f"{os.path.dirname(os.getcwd())}/output/cavity_models/best_model_path.txt",
            "r",
        )
        .read()
        .strip()
    )

    # Split into train and val
    dataloader_train_ds, dataloader_val_ds = train_val_split_ds(
        df_ddg, BATCH_SIZE_DS, DEVICE
    )

    # Define model
    cavity_model_net = CavityModel(get_latent=True).to(DEVICE)
    cavity_model_net.load_state_dict(
        torch.load(f"{best_cavity_model_path}", map_location=torch.device(DEVICE))
    )
    cavity_model_net.eval()
    loss_ds = torch.nn.L1Loss()

    # Train downstream model ensemble
    print("Starting downstream model ensemble training")
    for i in range(NUM_ENSEMBLE):
        # Initialize model with fixed seed
        model_idx = i
        ds_model_net = DownstreamModel().to(DEVICE)
        torch.manual_seed(seed=model_idx)
        torch.cuda.manual_seed(seed=model_idx)
        ds_model_net.apply(init_lin_weights)
        optimizer_ds = torch.optim.Adam(ds_model_net.parameters(), lr=LEARNING_RATE_DS)
        # Train model
        print(f"Training model: {model_idx+1}/{NUM_ENSEMBLE}")
        ds_train_val(
            df_ddg,
            dataloader_train_ds,
            dataloader_val_ds,
            cavity_model_net,
            ds_model_net,
            loss_ds,
            optimizer_ds,
            model_idx,
            EPOCHS_DS,
            DEVICE,
        )
        print(f"Finished training model: {model_idx+1}/{NUM_ENSEMBLE}")
    print("Finished downstream model ensemble training")

    # Plot learning curve with error bars
    learning_curve_ds_with_errorbars()

    # Load structure data
    pdb_filenames_test = {
        "Rosetta": sorted(
            glob.glob(
                f"{os.path.dirname(os.getcwd())}/data/test/Rosetta_10/structure/parsed/*coord*"
            )
        ),
        "Protein_G": sorted(
            glob.glob(
                f"{os.path.dirname(os.getcwd())}/data/test/Protein_G/structure/parsed/*coord*"
            )
        ),
        "ProTherm": sorted(
            glob.glob(
                f"{os.path.dirname(os.getcwd())}/data/test/ProTherm/targets/structure/parsed/*coord*"
            )
        ),
        "VAMP": sorted(
            glob.glob(
                f"{os.path.dirname(os.getcwd())}/data/test/VAMP/structure/parsed/*coord*"
            )
        ),
        "Homology": sorted(
            glob.glob(
                f"{os.path.dirname(os.getcwd())}/data/test/ProTherm/homology_models/structure/parsed/*coord*"
            )
        ),
        "Xstal_vs_AF2": sorted(
            glob.glob(
                f"{os.path.dirname(os.getcwd())}/data/test/Xstal_vs_AF2/structure/parsed/*coord*"
            )
        ),
    }

    # Load Rosetta ddGs
    rosetta_dict_test = {
        "Rosetta": pd.read_csv(
            f"{os.path.dirname(os.getcwd())}/data/test/Rosetta_10/ddG_Rosetta/ddg.csv"
        ),
        "Protein_G": pd.read_csv(
            f"{os.path.dirname(os.getcwd())}/data/test/Protein_G/ddG_Rosetta/ddg.csv"
        ),
        "ProTherm": pd.read_csv(
            f"{os.path.dirname(os.getcwd())}/data/test/ProTherm/targets/ddG_Rosetta/ddg.csv"
        ),
        "VAMP": pd.read_csv(
            f"{os.path.dirname(os.getcwd())}/data/test/VAMP/ddG_Rosetta/ddg.csv"
        ),
        "Homology": pd.read_csv(
            f"{os.path.dirname(os.getcwd())}/data/test/ProTherm/homology_models/ddG_Rosetta/ddg.csv"
        ),
        "Xstal_vs_AF2": pd.read_csv(
            f"{os.path.dirname(os.getcwd())}/data/test/Xstal_vs_AF2/ddG_Rosetta/ddg.csv"
        ),  # tmp
    }

    # Load experimental ddGs and VAMP-seq score data
    exp_dict_test = {
        "Protein_G": pd.read_csv(
            f"{os.path.dirname(os.getcwd())}/data/test/Protein_G/ddG_experimental/ddg.csv"
        ),
        "ProTherm": pd.read_csv(
            f"{os.path.dirname(os.getcwd())}/data/test/ProTherm/targets/ddG_experimental/ddg.csv"
        ),
        "VAMP": pd.read_csv(
            f"{os.path.dirname(os.getcwd())}/data/test/VAMP/VAMP/VAMP.csv"
        ),
        "Homology": pd.read_csv(
            f"{os.path.dirname(os.getcwd())}/data/test/ProTherm/homology_models/ddG_experimental/ddg.csv"
        ),
    }

    # Match structure, Rosetta and experimental data
    total_dict_test = {}

    for dataset_key, pdb_filenames in pdb_filenames_test.items():
        # Create temporary residue environment datasets to more easily match ddG data
        dataset_structure = ResidueEnvironmentsDataset(pdb_filenames, transformer=None)
        resenv_dataset = {}
        for resenv in dataset_structure:
            key = (
                f"{resenv.pdb_id}{resenv.chain_id}_{resenv.pdb_residue_number}"
                f"{index_to_one(resenv.restype_index)}"
            )
            resenv_dataset[key] = resenv

        # Populate Rosetta dataframes with wt ResidueEnvironment objects and wt and mt restype indices
        df_rosetta = rosetta_dict_test[dataset_key]
        df_rosetta = df_rosetta.rename(columns={"score": "score_rosetta"})
        n_rosetta_start = len(df_rosetta)
        populate_dfs_with_resenvs(df_rosetta, resenv_dataset)
        print(
            f"{n_rosetta_start-len(df_rosetta)} data points dropped when (inner) matching Rosetta and structure in: {dataset_key} data set."
        )

        # Add experimental ddGs to total dataframe
        if dataset_key in exp_dict_test.keys():
            df_exp = exp_dict_test[dataset_key]
            df_exp = df_exp.rename(columns={"score": "score_exp"})
            if dataset_key == "Homology":
                df_rosetta["pdbid_target"] = df_rosetta["pdbid"].str[:4]
                df_exp["pdbid_target"] = df_exp["pdbid"].str[:4]
                df_rosetta = df_rosetta.merge(
                    df_exp, on=["pdbid_target", "chainid", "variant"], how="inner"
                )
                df_rosetta["pdbid"] = df_rosetta["pdbid_x"]
                df_rosetta.drop(["pdbid_x", "pdbid_y"], axis=1)
            else:
                df_rosetta = df_rosetta.merge(
                    df_exp, on=["pdbid", "chainid", "variant"], how="inner"
                )
            print(
                f"{len(df_exp)-len(df_rosetta)} data points dropped when (inner) matching Rosetta and structural data with experimental data in: {dataset_key} data set."
            )
        df_total = df_rosetta

        # Remove disulfides
        if dataset_key == "Rosetta":
            n_df_start = len(df_total)
            df_total = remove_disulfides(df_total)
            print(
                f"{n_df_start-len(df_total)} data points removed due to disulfides in: {dataset_key} data set."
            )

        # Do Fermi transform
        df_total["score_rosetta_fermi"] = df_total["score_rosetta"].apply(
            fermi_transform
        )
        if "score_exp" in df_total:
            df_total["score_exp_fermi"] = df_total["score_exp"].apply(fermi_transform)

        # Add to total dict
        total_dict_test[dataset_key] = df_total

    # Add Protein G to ProTherm data set
    total_dict_test["ProTherm_with_Protein_G"] = pd.concat(
        [total_dict_test["ProTherm"], total_dict_test["Protein_G"]], axis=0
    )

    # Initialize models
    cavity_model_net = CavityModel(get_latent=True).to(DEVICE)
    best_cavity_model_path = (
        open(
            f"{os.path.dirname(os.getcwd())}/output/cavity_models/best_model_path.txt",
            "r",
        )
        .read()
        .strip()
    )
    cavity_model_net.load_state_dict(
        torch.load(f"{best_cavity_model_path}", map_location=torch.device(DEVICE))
    )
    cavity_model_net.eval()
    ds_model_net = DownstreamModel().to(DEVICE)

    # Make test set predictions
    for dataset_key, df_total in total_dict_test.items():
        print(f"Computing predictions for data set: {dataset_key}")

        # Compute predictions
        df_ml = ds_pred(
            cavity_model_net, ds_model_net, df_total, dataset_key, NUM_ENSEMBLE, DEVICE
        )

        # Merge and save data with predictions
        df_total = df_total.merge(
            df_ml, on=["pdbid", "chainid", "variant"], how="inner"
        )
        if dataset_key != "Homology":
            print(
                f"{len(df_ml)-len(df_total)} data points dropped when (inner) matching total data with ml predictions in: {dataset_key} data set."
            )
        df_total.to_csv(
            f"{os.path.dirname(os.getcwd())}/output/{dataset_key}/ml_preds.csv",
            index=False,
        )

        # Make visualizations
        if dataset_key == "Xstal_vs_AF2":
            compute_pdb_combo_corrs(df_total, dataset_key)
        else:
            scatter_plots(df_total, dataset_key)
            loss_analysis(df_total, dataset_key)
            if dataset_key == "Homology":
                homology_plot(df_total, dataset_key)
        print(f"Finished making predictions for data set: {dataset_key}")

    # Plot large cytosolic data set
    filenames = sorted(
        glob.glob(f"{os.path.dirname(os.getcwd())}/data/test/Cytosolic/*.csv")
    )
    df = pd.DataFrame()
    for file in filenames:
        df_file = pd.read_csv(file)
        df = pd.concat((df, df_file), axis=0)
    hist_plot_all(df, "Cytosolic")

    # Plot large cytosolic data set vs gnomAD and ClinVar
    plot_gnomad_clinvar()


if __name__ == "__main__":
    main()

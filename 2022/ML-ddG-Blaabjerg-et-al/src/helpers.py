import datetime
import glob
import itertools
import os
import pickle
import random
import sys
from typing import Dict, List, Union

import numpy as np
import pandas as pd
import pytz
import torch
from Bio.PDB.Polypeptide import index_to_one, one_to_index
from scipy.stats import pearsonr
from torch.nn.functional import softmax
from torch.utils.data import DataLoader, Dataset

from cavity_model import (
    CavityModel,
    DDGDataset,
    DDGToTensor,
    DownstreamModel,
    ResidueEnvironmentsDataset,
    ToTensor,
)
from PrismData import PrismParser, VariantData
from visualization import learning_curve_cavity, learning_curve_ds


def train_val_split_cavity(
    parsed_pdb_filenames,
    TRAIN_VAL_SPLIT,
    BATCH_SIZE,
    DEVICE,
):
    """
    Performs training and validation split of ResidueEnvironments.
    Note that we do the split on PDB level not on ResidueEnvironment level due to possible leakage.
    """
    n_train_pdbs = int(len(parsed_pdb_filenames) * TRAIN_VAL_SPLIT)
    filenames_train = parsed_pdb_filenames[:n_train_pdbs]
    filenames_val = parsed_pdb_filenames[n_train_pdbs:]
    to_tensor_transformer = ToTensor(DEVICE)

    dataset_train = ResidueEnvironmentsDataset(
        filenames_train, transformer=to_tensor_transformer
    )
    dataset_val = ResidueEnvironmentsDataset(
        filenames_val, transformer=to_tensor_transformer
    )

    dataloader_train = DataLoader(
        dataset_train,
        batch_size=BATCH_SIZE,
        shuffle=True,
        collate_fn=to_tensor_transformer.collate_cat,
        drop_last=True,
    )

    dataloader_val = DataLoader(
        dataset_val,
        batch_size=BATCH_SIZE,
        shuffle=False,
        collate_fn=to_tensor_transformer.collate_cat,
        drop_last=True,
    )

    print(
        f"Training data set includes {len(filenames_train)} pdbs with "
        f"{len(dataset_train)} environments."
    )
    print(
        f"Validation data set includes {len(filenames_val)} pdbs with "
        f"{len(dataset_val)} environments."
    )

    return dataloader_train, dataset_train, dataloader_val, dataset_val


def train_step(
    batch_x,
    batch_y,
    cavity_model_net,
    optimizer,
    loss_function,
):
    """
    Takes a training step with the cavity model
    """
    cavity_model_net.train()
    optimizer.zero_grad()
    batch_y_pred = cavity_model_net(batch_x)
    loss_batch = loss_function(batch_y_pred, torch.argmax(batch_y, dim=-1))
    loss_batch.backward()
    optimizer.step()
    return (batch_y_pred, loss_batch.detach().cpu().item())


def eval_loop(
    cavity_model_net,
    dataloader_val,
    loss_function,
):
    """
    Performs an eval loop with the cavity model
    """
    # Eval loop. Due to memory, we don't pass the whole eval set to the model
    labels_true_val = []
    labels_pred_val = []
    loss_batch_list_val = []
    cavity_model_net.eval()

    for batch_x_val, batch_y_val in dataloader_val:
        batch_y_pred_val = cavity_model_net(batch_x_val)

        loss_batch_val = loss_function(
            batch_y_pred_val, torch.argmax(batch_y_val, dim=-1)
        )
        loss_batch_list_val.append(loss_batch_val.detach().cpu().item())

        labels_true_val.append(torch.argmax(batch_y_val, dim=-1).detach().cpu().numpy())
        labels_pred_val.append(
            torch.argmax(batch_y_pred_val, dim=-1).detach().cpu().numpy()
        )

    acc_val = np.mean(
        (np.reshape(labels_true_val, -1) == np.reshape(labels_pred_val, -1))
    )
    loss_val = np.mean(loss_batch_list_val)

    return acc_val, loss_val


def train_loop(
    dataloader_train,
    dataloader_val,
    cavity_model_net,
    loss_function,
    optimizer,
    EPOCHS,
    PATIENCE_CUTOFF,
):
    current_best_epoch_idx = -1
    current_best_loss_val = 1e4
    patience = 0
    epoch_idx_to_model_path = {}

    acc_val_list = []
    acc_train_list = []
    loss_train_list = []
    """
    Trains the cavity model in a loop over EPOCHS
    """
    for epoch in range(EPOCHS):
        labels_true = []
        labels_pred = []
        loss_batch_list = []
        cavity_model_net.train()

        for batch_x, batch_y in dataloader_train:
            batch_y_pred, loss_batch = train_step(
                batch_x, batch_y, cavity_model_net, optimizer, loss_function
            )
            loss_batch_list.append(loss_batch)

            labels_true.append(torch.argmax(batch_y, dim=-1).detach().cpu().numpy())
            labels_pred.append(
                torch.argmax(batch_y_pred, dim=-1).detach().cpu().numpy()
            )

        # Train epoch metrics
        acc_train = np.mean(
            (np.reshape(labels_true, -1) == np.reshape(labels_pred, -1))
        )
        loss_train = np.mean(loss_batch_list)

        # Validation epoch metrics
        acc_val, loss_val = eval_loop(cavity_model_net, dataloader_val, loss_function)
        acc_val_list.append(acc_val)
        acc_train_list.append(acc_train)
        loss_train_list.append(loss_train)

        print(
            f"Epoch {epoch:2d}. Train loss: {loss_train:5.3f}. "
            f"Train Acc: {acc_train:4.2f}. Val loss: {loss_val:5.3f}. "
            f"Val Acc {acc_val:4.2f}"
        )

        # Save model
        model_path = f"{os.path.dirname(os.getcwd())}/output/cavity_models/cavity_model_{epoch:02d}.pt"
        epoch_idx_to_model_path[epoch] = model_path
        torch.save(cavity_model_net.state_dict(), model_path)

        # Early stopping
        if loss_val < current_best_loss_val:
            current_best_loss_val = loss_val
            current_best_epoch_idx = epoch
            patience = 0
        else:
            patience += 1
        if patience > PATIENCE_CUTOFF:
            print(f"Early stopping activated.")
            break

    learning_curve_cavity(acc_val_list, acc_train_list, loss_train_list)
    best_model_path = epoch_idx_to_model_path[current_best_epoch_idx]
    print(
        f"Best epoch idx: {current_best_epoch_idx} with validation loss: "
        f"{current_best_loss_val:5.3f} and model_path: "
        f"{best_model_path}"
    )
    return best_model_path


def populate_dfs_with_resenvs(ddg_data, resenv_dataset):
    """
    Populates ddG dfs with the wt ResidueEnvironment objects.
    """
    print(
        "Dropping data points where residue is not defined in structure "
        f"or due to missing parsed pdb file"
    )
    # Add wt residue environments to standard ddg data dataframes

    resenvs_ddg_data = []
    for idx, row in ddg_data.iterrows():
        resenv_key = (
            f"{row['pdbid']}{row['chainid']}_"
            f"{row['variant'][1:-1]}{row['variant'][0]}"
        )
        try:
            resenv = resenv_dataset[resenv_key]
            resenvs_ddg_data.append(resenv)
        except KeyError:
            print(resenv_key + " " "Could not be found in structure data")
            resenvs_ddg_data.append(np.nan)
    ddg_data["resenv"] = resenvs_ddg_data
    n_datapoints_before = ddg_data.shape[0]
    ddg_data.dropna(inplace=True)
    n_datapoints_after = ddg_data.shape[0]
    print(
        f"dropped {n_datapoints_before - n_datapoints_after:4d} / "
        f"{n_datapoints_before:4d} data points from dataset."
    )

    # Load PDB amino acid frequencies used to approximate unfolded states
    pdb_nlfs = -np.log(
        np.load(
            f"{os.path.dirname(os.getcwd())}/data/train/cavity/pdb_frequencies.npz"
        )["frequencies"]
    )

    # Add wt and mt idxs and freqs to df
    ddg_data["wt_idx"] = ddg_data.apply(
        lambda row: one_to_index(row["variant"][0]), axis=1
    )
    ddg_data["mt_idx"] = ddg_data.apply(
        lambda row: one_to_index(row["variant"][-1]), axis=1
    )
    ddg_data["wt_nlf"] = ddg_data.apply(lambda row: pdb_nlfs[row["wt_idx"]], axis=1)
    ddg_data["mt_nlf"] = ddg_data.apply(lambda row: pdb_nlfs[row["mt_idx"]], axis=1)


def remove_disulfides(ddg_data):
    """
    Removes residues known to be involved in disulfide bridges from the data
    """
    disulfide_dict = {
        "6B0N": [
            "C23",
            "C42",
            "C87",
            "C94",
            "C99",
            "C116",
            "C154",
            "C163",
            "C176",
            "C186",
            "C197",
            "C205",
            "C254",
            "C288",
            "C335",
            "C342",
            "C365",
            "C392",
            "C544",
            "C550",
        ],
        "1ZG4": ["C52", "C98"],
        "2YP7": [
            "C7",
            "C45",
            "C57",
            "C69",
            "C90",
            "C132",
            "C270",
            "C274",
            "C298",
            "C453",
            "C460",
            "C464",
        ],
        "1PS1": ["C126", "C134"],
        "6M17": ["C113", "C121", "C341", "C510", "C522"],
        "2R7E": ["C322", "C241", "C179", "C153", "C598", "C679", "C496", "C522"],
        "4WVP": ["C179", "C122", "C158", "C152", "C42", "C26", "C194", "C169"],
    }
    for pdb in disulfide_dict:
        for res in disulfide_dict[pdb]:
            ddg_data = ddg_data[
                ~((ddg_data["pdbid"] == pdb) & (ddg_data["variant"].str.contains(res)))
            ]

    return ddg_data


def fermi_transform(x):
    """
    Applies a Fermi transformation to values
    """
    alpha = 3.0
    beta = 0.4
    y = 1.0 / (1.0 + np.exp(-beta * (x - alpha)))
    return y


def inverse_fermi_transform(x):
    """
    Applies an inverse Fermi transformation to values
    """
    alpha = 3.0
    beta = 0.4
    EPS = 10.0 ** (-12)
    y = 0.0
    if x == 1.0:
        y = 40.0
    elif x > 0.0 and 1.0 > x:
        y = (alpha * beta - np.log(-1.0 + 1.0 / x + EPS)) / beta
    elif x == 0.0:
        y = -40.0
    return y


def get_ddg_dataloader(ddg_data, data_type, BATCH_SIZE, DEVICE):
    """
    Returns a ddG dataloader based on input data and data type
    """
    # Define ddG data set
    ddg_dataset = DDGDataset(ddg_data, transformer=DDGToTensor(data_type, DEVICE))

    # Define case-dependent data loader parameters
    shuffle = False
    if data_type == "train":
        BATCH_SIZE = BATCH_SIZE
        shuffle = True
        drop_last = True
    elif data_type == "val":
        BATCH_SIZE = BATCH_SIZE
        shuffle = False
        drop_last = True
    elif data_type == "test" or data_type == "pred":
        BATCH_SIZE = 100  # Worst-case max gpu memory on small machine
        shuffle = False
        drop_last = False

    # Define data loader
    ddg_dataloader = DataLoader(
        ddg_dataset,
        batch_size=BATCH_SIZE,
        shuffle=shuffle,
        drop_last=drop_last,
        collate_fn=DDGToTensor(data_type, DEVICE).collate_multi,
    )

    return ddg_dataloader


def train_val_split_ds(ddg_data, BATCH_SIZE, DEVICE):
    """
    Performs training and validation split of ddG dataframe.
    Note that we do the split on PDB level not on ResidueEnvironment level due to possible leakage.
    """
    parsed_pdb_filenames = ddg_data["pdbid"].unique()
    filenames_val = np.array(
        ["2GRN", "3UFJ", "6R5K", "6B0N", "4QO1", "1BRW", "2OCP", "1DO6", "1ZG4", "1ARZ"]
    )
    filenames_train = np.setdiff1d(parsed_pdb_filenames, filenames_val)

    ddg_data_train = ddg_data[ddg_data["pdbid"].isin(filenames_train)]
    ddg_data_val = ddg_data[ddg_data["pdbid"].isin(filenames_val)]

    dataloader_train = get_ddg_dataloader(ddg_data_train, "train", BATCH_SIZE, DEVICE)
    dataloader_val = get_ddg_dataloader(ddg_data_val, "val", BATCH_SIZE, DEVICE)

    print(
        f"Training data set includes {len(filenames_train)} pdbs with "
        f"{len(ddg_data_train)} mutations."
    )
    print(f"Training PDBs are: {filenames_train}.")
    print(
        f"Validation data set includes {len(filenames_val)} pdbs with "
        f"{len(ddg_data_val)} mutations."
    )
    print(f"Validation PDBs are: {filenames_val}.")

    return dataloader_train, dataloader_val


def init_lin_weights(m):
    if type(m) == torch.nn.Linear:
        torch.nn.init.kaiming_uniform_(m.weight)
        m.bias.data.fill_(0.01)


def ds_train_val(
    df_ddg,
    dataloader_train,
    dataloader_val,
    cavity_model_net,
    ds_model_net,
    loss_func,
    optimizer,
    model_idx,
    EPOCHS,
    DEVICE,
):
    """
    Trains the downstream model in a loop over EPOCHS.
    """

    # Set seed
    np.random.seed(model_idx)
    random.seed(model_idx)
    torch.manual_seed(model_idx)
    torch.cuda.manual_seed(model_idx)
    torch.cuda.manual_seed_all(model_idx)

    # Initialize
    train_loss_list = []
    val_loss_list = []
    pearson_r_list = []

    for epoch in range(EPOCHS):

        print(f"Epoch: {epoch+1}/{EPOCHS}")

        # Initialize
        train_loss_batch_list = torch.empty(0, dtype=torch.float32).to(DEVICE)
        val_loss_batch_list = torch.empty(0, dtype=torch.float32).to(DEVICE)
        val_ddg_fermi = torch.empty(0, dtype=torch.float32).to(DEVICE)
        val_ddg_fermi_pred = torch.empty(0, dtype=torch.float32).to(DEVICE)
        val_ddg = torch.empty(0, dtype=torch.float32).to(DEVICE)

        # Train loop
        ds_model_net.train()

        for _, _, _, x_cavity_batch, x_ds_batch, ddg_fermi_batch, _ in dataloader_train:

            # Initialize optimizer
            optimizer.zero_grad()

            # Compute predictions
            cavity_pred_batch = cavity_model_net(x_cavity_batch)
            ddg_fermi_batch_pred = ds_model_net(
                torch.cat((cavity_pred_batch, x_ds_batch), 1)
            )
            loss_batch = loss_func(ddg_fermi_batch_pred, ddg_fermi_batch)
            loss_batch.backward()
            optimizer.step()

            # Append to epoch
            train_loss_batch_list = torch.cat(
                (train_loss_batch_list, loss_batch.detach().reshape(-1))
            )

        # Val loop
        ds_model_net.eval()
        with torch.no_grad():

            for (
                _,
                _,
                _,
                val_x_cavity_batch,
                val_x_ds_batch,
                val_ddg_fermi_batch,
                val_ddg_batch,
            ) in dataloader_val:

                # Compute predictions
                val_cavity_pred_batch = cavity_model_net(val_x_cavity_batch)
                val_ddg_fermi_pred_batch = ds_model_net(
                    torch.cat((val_cavity_pred_batch, val_x_ds_batch), 1)
                )
                val_loss_batch = loss_func(
                    val_ddg_fermi_pred_batch, val_ddg_fermi_batch
                )

                # Append to epoch
                val_loss_batch_list = torch.cat(
                    (val_loss_batch_list, val_loss_batch.reshape(-1))
                )
                val_ddg_fermi = torch.cat(
                    (val_ddg_fermi, val_ddg_fermi_batch.reshape(-1)), 0
                )
                val_ddg_fermi_pred = torch.cat(
                    (val_ddg_fermi_pred, val_ddg_fermi_pred_batch.reshape(-1)), 0
                )
                val_ddg = torch.cat((val_ddg, val_ddg_batch.reshape(-1)), 0)

        # Compute epoch metrics
        train_loss_list.append(train_loss_batch_list.mean().cpu().item())
        val_loss_list.append(val_loss_batch_list.mean().cpu().item())
        val_ddg = val_ddg.detach().cpu()
        val_ddg_fermi_pred = val_ddg_fermi_pred.detach().cpu()
        val_ddg_pred = val_ddg_fermi_pred.apply_(lambda x: inverse_fermi_transform(x))
        pearson_r_list.append(pearsonr(val_ddg_pred.numpy(), val_ddg.numpy())[0])

    # Print learning curve
    learning_curve_ds(pearson_r_list, val_loss_list, train_loss_list, model_idx)

    # Save results and model
    lc_results = (pearson_r_list, val_loss_list, train_loss_list, model_idx)
    with open(
        f"{os.path.dirname(os.getcwd())}/output/ds_models/ds_model_{model_idx}/lc_results.pkl",
        "wb",
    ) as f:
        pickle.dump(lc_results, f)
    torch.save(
        ds_model_net.state_dict(),
        f"{os.path.dirname(os.getcwd())}/output/ds_models/ds_model_{model_idx}/model.pt",
    )


def ds_pred(
    cavity_model_net, ds_model_net, df_total, dataset_key, NUM_ENSEMBLE, DEVICE
):
    """
    Makes predictions using the cavity and downstream model (RaSP model)
    """

    # Make data loader
    dataloader = get_ddg_dataloader(df_total, "pred", None, DEVICE)

    # Make predictions
    pdbid = []
    chainid = []
    variant = []
    ddg_fermi_pred = torch.empty(0, 1, dtype=torch.float32).to(DEVICE)

    with torch.no_grad():
        for (
            pdbid_batch,
            chainid_batch,
            variant_batch,
            x_cavity_batch,
            x_ds_batch,
        ) in dataloader:

            # Initialize
            ddg_fermi_pred_batch_ensemble = torch.empty(len(variant_batch), 0).to(
                DEVICE
            )

            # Load ds models in a load
            for i in range(NUM_ENSEMBLE):
                model_idx = i
                ds_model_net.load_state_dict(
                    torch.load(
                        f"{os.path.dirname(os.getcwd())}/output/ds_models/ds_model_{model_idx}/model.pt"
                    )
                )
                ds_model_net.eval()

                # Compute predictions
                cavity_pred_batch = cavity_model_net(x_cavity_batch)
                ddg_fermi_pred_batch = ds_model_net(
                    torch.cat((cavity_pred_batch, x_ds_batch), 1)
                )
                ddg_fermi_pred_batch_ensemble = torch.cat(
                    (ddg_fermi_pred_batch_ensemble, ddg_fermi_pred_batch), 1
                )

            # Take median of ensemble predictions
            ddg_fermi_pred_batch = torch.median(
                ddg_fermi_pred_batch_ensemble, 1, keepdim=True
            )[0]

            # Append to epoch results
            pdbid += pdbid_batch
            chainid += chainid_batch
            variant += variant_batch
            ddg_fermi_pred = torch.cat((ddg_fermi_pred, ddg_fermi_pred_batch), 0)

    # Repack data to df
    df_ml = pd.DataFrame(ddg_fermi_pred.cpu().numpy(), columns=["score_ml_fermi"])
    df_ml["score_ml"] = df_ml["score_ml_fermi"].apply(
        lambda x: inverse_fermi_transform(x)
    )
    df_ml.insert(loc=0, column="pdbid", value=np.array(pdbid))
    df_ml.insert(loc=1, column="chainid", value=np.array(chainid))
    df_ml.insert(loc=2, column="variant", value=np.array(variant))
    return df_ml


def cavity_to_prism(df, pdbid, chainid, seq, prism_file):
    """
    Saves a ddG dataframe to a PRISM-format csv file
    """

    # Get relevant columns
    df_prism = df[["variant", "score_ml_fermi", "score_ml"]]

    # Get current CPH time
    timestamp = datetime.datetime.now(pytz.timezone("Europe/Copenhagen")).strftime(
        "%Y-%m-%d %H:%M"
    )

    metadata = {
        "version": 1,
        "protein": {
            "name": "Unknown",
            "organism": "Unknown",
            "uniprot": "Unknown",
            "sequence": seq,
            "pdb": pdbid,
            "chain": chainid,
        },
        "cavity": {
            "version": 1,
        },
        "columns": {
            "score_ml_fermi": "Normalized cavity model ddG predictions (range [0;1])",
            "score_ml": "Cavity model ddG predictions",
        },
        "created": {f"{timestamp} (CPH time) - lasse.blaabjerg@bio.ku.dk"},
    }

    # Write data
    dataset = VariantData(metadata, df_prism)
    comment = ""
    parser = PrismParser()
    parser.write(prism_file, dataset, comment_lines=comment)


def get_seq_from_variant(df):
    """
    Fetches the wt amino acid sequence from a ddG dataframe
    """

    df["pos"] = df["variant"].str.extract("(\d+)").astype(int)
    df["wt"] = df["variant"].str.extract("(^[a-zA-Z]+)")
    seq = ""
    for i in range(1, max(df["pos"] + 1)):
        if i in df["pos"].unique():
            seq += df["wt"][df["pos"] == i].iloc[0]
        else:
            seq += "X"
    return seq


def compute_pdb_combo_corrs(df, dataset_key):
    """
    Computes correlation coefficients between ddGs from a list of PDBs
    """

    parser = PrismParser()

    df["pos"] = df["variant"].str.extract("(\d+)").astype(int)

    combos = [
        ("6V0P", "4GQB"),
        ("6V0P", "AFAR"),
        ("4GQB", "AFAR"),
        ("6B6U", "6NU5"),
        ("6B6U", "AFPK"),
        ("6NU5", "AFPK"),
        ("4Y08", "4OYN"),
        ("4Y08", "AFFH"),
        ("4OYN", "AFFH"),
    ]

    results = []

    for combo in combos:
        df_pdb_1 = df[df["pdbid"] == combo[0]].reset_index()
        df_pdb_2 = df[df["pdbid"] == combo[1]].reset_index()

        pdbid_1 = df_pdb_1["pdbid"][0]
        chainid_1 = df_pdb_1["chainid"][0]
        seq_1 = get_seq_from_variant(df_pdb_1)
        filename_1 = f"{os.path.dirname(os.getcwd())}/output/{dataset_key}/prism_cavity_pred_{pdbid_1}"
        cavity_to_prism(df_pdb_1, pdbid_1, chainid_1, seq_1, filename_1)
        df_pdb_1_prism = parser.read(filename_1)

        pdbid_2 = df_pdb_2["pdbid"][0]
        chainid_2 = df_pdb_2["chainid"][0]
        seq_2 = get_seq_from_variant(df_pdb_2)
        filename_2 = f"{os.path.dirname(os.getcwd())}/output/{dataset_key}/prism_cavity_pred_{pdbid_2}"
        cavity_to_prism(df_pdb_2, pdbid_2, chainid_2, seq_2, filename_2)
        df_pdb_2_prism = parser.read(filename_2)

        # Merge
        df_total = df_pdb_1_prism.merge([df_pdb_2_prism]).dataframe.dropna()
        x = df_total["score_ml_00"]
        y = df_total["score_ml_01"]

        y = y[((x >= -1) & (x <= 7)) & ((y >= -1) & (y <= 7))]
        x = x[((x >= -1) & (x <= 7)) & ((y >= -1) & (y <= 7))]

        pearson_r = pearsonr(x, y)[0]

        print(
            f"Pearson r of RaSP predictions between {combo[0]} and {combo[1]} is: {pearson_r}"
        )
        results.append((combo[0], combo[1], pearson_r))

    with open(
        f"{os.path.dirname(os.getcwd())}/output/{dataset_key}/combo_corr_results.pkl",
        "wb",
    ) as f:
        pickle.dump(results, f)

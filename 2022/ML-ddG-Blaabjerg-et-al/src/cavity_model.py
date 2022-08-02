import glob
import os
import random
from typing import Callable, List, Union

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset

__all__ = [
    "ResidueEnvironment",
    "ResidueEnvironmentsDataset",
    "ToTensor",
    "CavityModel",
    "DownstreamModel",
    "DDGDataset",
    "DDGToTensor",
]


class ResidueEnvironment:
    """
    Residue environment class used to hold necessarry information about the
    atoms of the environment such as atomic coordinates, atom types and the
    class of the missing central amino acid.

    Parameters
    ----------
    xyz_coords: np.ndarray
        Numpy array with shape (n_atoms, 3) containing the x, y, z coordinates.
    atom_types: np.ndarray
        1D numpy array containing the atom types. Integer values in range(6).
    restypes_onehot: np.ndarray
        Numpy array with shape (n_atoms, 21) containing the amino acid
        class of the missing amino acid
    chain_id: str
        Chain id associated to ResidueEnvironment object
    pdb_residue_number: int
        Residue number associated with the ResidueEnvironment object
    pdb_id: str
        PDBID associated with the ResidueEnvironment object
    """

    def __init__(
        self,
        xyz_coords,
        atom_types,
        restype_onehot,
        chain_id,
        pdb_residue_number,
        pdb_id,
    ):
        self.xyz_coords = xyz_coords
        self.atom_types = atom_types
        self.restype_onehot = restype_onehot
        self.restype_index = np.argmax(self.restype_onehot)
        self.chain_id = chain_id
        self.pdb_residue_number = pdb_residue_number
        self.pdb_id = pdb_id

    def __repr__(self):
        return (
            f"<ResidueEnvironment with {self.xyz_coords.shape[0]} atoms. "
            f"pdb_id: {self.pdb_id}, "
            f"chain_id: {self.chain_id}, "
            f"pdb_residue_number: {self.pdb_residue_number}, "
            f"restype_index: {self.restype_index}>"
        )


class ResidueEnvironmentsDataset(Dataset):
    """
    Residue environment dataset class

    Parameters
    ----------
    input_data: Union[List[str], List[ResidueEnvironment]]
        List of parsed pdb filenames in .npz format or list of
        ResidueEnvironment objects
    transform: Callable
        A to-tensor transformer class
    """

    def __init__(
        self,
        input_data,
        transformer,
    ):

        if all(isinstance(x, ResidueEnvironment) for x in input_data):
            self.res_env_objects = input_data
        elif all(isinstance(x, str) for x in input_data):
            self.res_env_objects = self.parse_envs(input_data)
        else:
            raise ValueError(
                "Input data is not of type" "Union[List[str], List[ResidueEnvironment]]"
            )

        self.transformer = transformer

    def __len__(self):
        return len(self.res_env_objects)

    def __getitem__(self, idx):
        sample = self.res_env_objects[idx]
        if self.transformer:
            sample = self.transformer(sample)
        return sample

    def parse_envs(self, npz_filenames):
        res_env_objects = []
        for i in range(len(npz_filenames)):
            coordinate_features = np.load(npz_filenames[i])
            atom_coords_prot_seq = coordinate_features["positions"]
            restypes_onehots_prot_seq = coordinate_features["aa_onehot"]
            selector_prot_seq = coordinate_features["selector"]
            atom_types_flattened = coordinate_features["atom_types_numeric"]

            chain_ids = coordinate_features["chain_ids"]
            pdb_residue_numbers = coordinate_features["residue_numbers"]
            chain_boundary_indices = coordinate_features["chain_boundary_indices"]
            N_residues = selector_prot_seq.shape[0]
            if "bin" in npz_filenames[i]:  # Check if file belongs to Homology data set
                pdb_id = os.path.basename(npz_filenames[i])[:19]
            else:
                pdb_id = os.path.basename(npz_filenames[i])[:4]

            for resi_i in range(N_residues):
                selector = selector_prot_seq[resi_i]
                selector_masked = selector[selector > -1]  # Remove Filler
                coords_mask = (
                    atom_coords_prot_seq[resi_i, :, 0] != -99.0
                )  # Remove filler
                coords = atom_coords_prot_seq[resi_i][coords_mask]
                atom_types = atom_types_flattened[selector_masked]
                restype_onehot = restypes_onehots_prot_seq[resi_i]

                pdb_residue_number = int(pdb_residue_numbers[resi_i])
                # Locate chain id
                for j in range(len(chain_ids)):
                    chain_boundary_0 = chain_boundary_indices[j]
                    chain_boundary_1 = chain_boundary_indices[j + 1]
                    if resi_i in range(chain_boundary_0, chain_boundary_1):
                        chain_id = str(chain_ids[j])
                        break

                res_env_objects.append(
                    ResidueEnvironment(
                        coords,
                        atom_types,
                        restype_onehot,
                        chain_id,
                        pdb_residue_number,
                        pdb_id,
                    )
                )

        return res_env_objects


class ToTensor:
    """
    To-tensor transformer

    Parameters
    ----------
    DEVICE: str
        Either "cuda" (gpu) or "cpu". Is set-able.
    """

    def __init__(self, DEVICE):
        self.device = DEVICE

    def __call__(self, sample):
        """
        Converts single ResidueEnvironment object into x_ and y_
        """

        sample_env = np.hstack(
            [np.reshape(sample.atom_types, [-1, 1]), sample.xyz_coords]
        )

        return {
            "x_": torch.tensor(sample_env, dtype=torch.float32).to(self.device),
            "y_": torch.tensor(np.array(sample.restype_onehot), dtype=torch.float32).to(
                self.device
            ),
        }

    def collate_cat(self, batch):
        """
        Collate method used by the dataloader to collate a
        batch of ResidueEnvironment objects.
        """
        target = torch.cat([torch.unsqueeze(b["y_"], 0) for b in batch], dim=0)

        # To collate the input, we need to add a column which
        # specifies the environtment each atom belongs to
        env_id_batch = []
        for i, b in enumerate(batch):
            n_atoms = b["x_"].shape[0]
            env_id_arr = torch.zeros(n_atoms, dtype=torch.float32).to(self.device) + i
            env_id_batch.append(
                torch.cat([torch.unsqueeze(env_id_arr, 1), b["x_"]], dim=1)
            )
        data = torch.cat(env_id_batch, dim=0)

        return data, target


class CavityModel(torch.nn.Module):
    """
    3D convolutional neural network to missing amino acid classification

    Parameters
    ----------
    DEVICE: str
        Either "cuda" (gpu) or "cpu". Is set-able.
    n_atom_types: int
        Number of atom types. (C, H, N, O, S, P)
    bins_per_angstrom: float
        Number of grid points per Anstrom.
    grid_dim: int
        Grid dimension
    sigma: float
        Standard deviation used for gaussian blurring
    """

    def __init__(
        self,
        get_latent=False,
        n_atom_types=6,
        bins_per_angstrom=1.0,
        grid_dim=18,
        sigma=0.6,
    ):
        super().__init__()

        self.n_atom_types = n_atom_types
        self.get_latent = get_latent
        self.bins_per_angstrom = bins_per_angstrom
        self.grid_dim = grid_dim
        self.sigma = sigma
        self.sigma_p = self.sigma * self.bins_per_angstrom
        self.model()

    def lin_spacing(self):
        lin_spacing = np.linspace(
            start=-self.grid_dim / 2 * self.bins_per_angstrom
            + self.bins_per_angstrom / 2,
            stop=self.grid_dim / 2 * self.bins_per_angstrom
            - self.bins_per_angstrom / 2,
            num=self.grid_dim,
        )
        return lin_spacing

    def model(self):
        self.xx, self.yy, self.zz = torch.tensor(
            np.meshgrid(
                self.lin_spacing(),
                self.lin_spacing(),
                self.lin_spacing(),
                indexing="ij",
            ),
            dtype=torch.float32,
        )

        self.conv1 = torch.nn.Sequential(
            torch.nn.Conv3d(6, 16, kernel_size=(3, 3, 3), stride=1, padding=1),
            torch.nn.MaxPool3d(
                kernel_size=(2, 2, 2), stride=(2, 2, 2), padding=(0, 0, 0)
            ),
            torch.nn.BatchNorm3d(16),
            torch.nn.LeakyReLU(),
        )
        self.conv2 = torch.nn.Sequential(
            torch.nn.Conv3d(16, 32, kernel_size=(3, 3, 3), stride=1, padding=0),
            torch.nn.MaxPool3d(
                kernel_size=(2, 2, 2), stride=(2, 2, 2), padding=(0, 0, 0)
            ),
            torch.nn.BatchNorm3d(32),
            torch.nn.LeakyReLU(),
        )
        self.conv3 = torch.nn.Sequential(
            torch.nn.Conv3d(32, 64, kernel_size=(3, 3, 3), stride=1, padding=1),
            torch.nn.MaxPool3d(
                kernel_size=(2, 2, 2), stride=(2, 2, 2), padding=(0, 0, 0)
            ),
            torch.nn.BatchNorm3d(64),
            torch.nn.LeakyReLU(),
            torch.nn.Flatten(),
        )
        self.dense1 = torch.nn.Sequential(
            torch.nn.Linear(in_features=64, out_features=128),
            torch.nn.BatchNorm1d(128),
            torch.nn.LeakyReLU(),
        )
        self.dense2 = torch.nn.Linear(in_features=128, out_features=100)
        self.dense3 = torch.nn.Linear(in_features=100, out_features=20)

    def forward(self, x):
        x = self.gaussian_blurring(x)
        x = self.conv1(x)
        x = self.conv2(x)
        x = self.conv3(x)
        x = self.dense1(x)
        x = self.dense2(x)
        if self.get_latent == True:
            return x
        elif self.get_latent == False:
            x = torch.nn.functional.leaky_relu(x)
            x = self.dense3(x)
            return x

    def gaussian_blurring(self, x):
        """
        Method that takes 2d torch.Tensor describing the atoms of the batch.

        Parameters
        ----------
        x: torch.Tensor
            Tensor for shape (n_atoms, 5). Each row represents an atom, where:
                column 0 describes the environment of the batch the
                atom belongs to
                column 1 describes the atom type
                column 2,3,4 are the x, y, z coordinates, respectively

        Returns
        -------
        fields_torch: torch.Tensor
            Represents the structural environment with gaussian blurring
            and has shape (-1, self.grid_dim, self.grid_dim, self.grid_dim).
        """
        current_batch_size = torch.unique(x[:, 0]).shape[0]
        fields_torch = torch.zeros(
            (
                current_batch_size,
                self.n_atom_types,
                self.grid_dim,
                self.grid_dim,
                self.grid_dim,
            )
        ).to(x.device)
        for j in range(self.n_atom_types):
            mask_j = x[:, 1] == j
            atom_type_j_data = x[mask_j]
            if atom_type_j_data.shape[0] > 0:
                pos = atom_type_j_data[:, 2:]
                density = torch.exp(
                    -(
                        (torch.reshape(self.xx.to(x.device), [-1, 1]) - pos[:, 0]) ** 2
                        + (torch.reshape(self.yy.to(x.device), [-1, 1]) - pos[:, 1])
                        ** 2
                        + (torch.reshape(self.zz.to(x.device), [-1, 1]) - pos[:, 2])
                        ** 2
                    )
                    / (2 * self.sigma_p ** 2)
                )

                # Normalize each atom to 1
                density /= torch.sum(density, dim=0)

                # Since column 0 of atom_type_j_data is sorted
                # I can use a trick to detect the boundaries based
                # on the change from one value to another.
                change_mask_j = (
                    atom_type_j_data[:, 0][:-1] != atom_type_j_data[:, 0][1:]
                )

                # Add begin and end indices
                ranges_i = torch.cat(
                    [
                        torch.tensor([0]),
                        torch.arange(atom_type_j_data.shape[0] - 1)[change_mask_j] + 1,
                        torch.tensor([atom_type_j_data.shape[0]]),
                    ]
                )

                # Fill tensor
                for i in range(ranges_i.shape[0]):
                    if i < ranges_i.shape[0] - 1:
                        index_0, index_1 = ranges_i[i], ranges_i[i + 1]
                        fields = torch.reshape(
                            torch.sum(density[:, index_0:index_1], dim=1),
                            [self.grid_dim, self.grid_dim, self.grid_dim],
                        )
                        fields_torch[i, j, :, :, :] = fields
        return fields_torch


class DownstreamModel(torch.nn.Module):
    """
    Downstream FC neural network.
    """

    def __init__(self):
        super().__init__()

        # Model
        self.lin1 = torch.nn.Sequential(
            torch.nn.Linear(142, 128),
            torch.nn.BatchNorm1d(128),
            torch.nn.LeakyReLU(),
        )
        self.lin2 = torch.nn.Sequential(
            torch.nn.Linear(128, 64),
            torch.nn.BatchNorm1d(64),
            torch.nn.LeakyReLU(),
        )
        self.lin3 = torch.nn.Sequential(
            torch.nn.Linear(64, 16),
            torch.nn.BatchNorm1d(16),
            torch.nn.LeakyReLU(),
        )
        self.lin4 = torch.nn.Linear(in_features=16, out_features=1)
        self.sigmoid = torch.nn.Sigmoid()

    # Forward
    def forward(self, x):
        x = self.lin1(x)
        x = self.lin2(x)
        x = self.lin3(x)
        x = self.lin4(x)
        x = self.sigmoid(x)
        return x


class DDGDataset(Dataset):
    """
    ddG dataset
    """

    def __init__(
        self,
        df,
        transformer=None,
    ):

        self.df = df
        self.transformer = transformer

    def __len__(self):
        return self.df.shape[0]

    def __getitem__(self, idx):
        sample = self.df.iloc[idx]
        if self.transformer:
            sample = self.transformer(sample)
        return sample


class DDGToTensor:
    """
    To-tensor transformer for ddG dataframe data

    Parameters
    ----------
    DEVICE: str
        Either "cuda" (gpu) or "cpu". Is set-able.
    """

    def __init__(self, flag, DEVICE):
        self.device = DEVICE
        self.flag = flag

    def __call__(self, sample):
        # Cavity
        atom_types = sample["resenv"].atom_types
        xyz_coords = sample["resenv"].xyz_coords

        x_cavity = np.hstack([np.reshape(atom_types, [-1, 1]), xyz_coords])

        # Downstream
        wt_onehot = np.zeros(20)
        wt_onehot[sample["wt_idx"]] = 1.0
        mt_onehot = np.zeros(20)
        mt_onehot[sample["mt_idx"]] = 1.0
        wt_freq = [sample["wt_nlf"]]
        mt_freq = [sample["mt_nlf"]]
        x_ds = np.hstack([wt_onehot, mt_onehot, wt_freq, mt_freq])
        pdbid = [sample["pdbid"]]
        chainid = [sample["chainid"]]
        variant = [sample["variant"]]

        if self.flag != "pred":
            ddg_fermi = [sample["score_fermi"]]
            ddg = [sample["score"]]
            return {
                "pdbid": pdbid,
                "chainid": chainid,
                "variant": variant,
                "x_cavity": torch.tensor(x_cavity, dtype=torch.float32).to(self.device),
                "x_ds": torch.tensor(x_ds, dtype=torch.float32).to(self.device),
                "score_fermi": torch.tensor(ddg_fermi, dtype=torch.float32).to(
                    self.device
                ),
                "score": torch.tensor(ddg, dtype=torch.float32).to(self.device),
            }
        else:
            return {
                "pdbid": pdbid,
                "chainid": chainid,
                "variant": variant,
                "x_cavity": torch.tensor(x_cavity, dtype=torch.float32).to(self.device),
                "x_ds": torch.tensor(x_ds, dtype=torch.float32).to(self.device),
            }

    def collate_multi(self, batch):
        """
        Collate method used by the dataloader to collate a
        batch of multi data.
        """
        # To collate the input, we need to add a column which
        # specifies the environment each atom belongs to
        env_id_batch = []
        for i, b in enumerate(batch):
            n_atoms = b["x_cavity"].shape[0]
            env_id_arr = torch.zeros(n_atoms, dtype=torch.float32) + i
            env_id_arr = env_id_arr.to(self.device)
            env_id_batch.append(
                torch.cat([torch.unsqueeze(env_id_arr, 1), b["x_cavity"]], dim=1)
            )
        x_cavity_batch = torch.cat(env_id_batch, dim=0)

        # Downstream
        pdbid_batch = [b["pdbid"] for b in batch]
        chainid_batch = [b["chainid"] for b in batch]
        variant_batch = [b["variant"] for b in batch]
        x_ds_batch = torch.cat([torch.unsqueeze(b["x_ds"], 0) for b in batch], dim=0)

        if self.flag != "pred":
            ddg_fermi_batch = torch.cat(
                [torch.unsqueeze(b["score_fermi"], 0) for b in batch], dim=0
            )
            ddg_batch = torch.cat(
                [torch.unsqueeze(b["score"], 0) for b in batch], dim=0
            )
            return (
                pdbid_batch,
                chainid_batch,
                variant_batch,
                x_cavity_batch,
                x_ds_batch,
                ddg_fermi_batch,
                ddg_batch,
            )
        else:
            return pdbid_batch, chainid_batch, variant_batch, x_cavity_batch, x_ds_batch

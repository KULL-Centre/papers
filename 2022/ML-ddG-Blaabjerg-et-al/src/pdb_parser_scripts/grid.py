import enum

import Bio
import numpy as np

# Enum for different typs of coordinate system
# DO NOT changes the numbers, as the number value is stored as a feature in files
CoordinateSystem = enum.Enum(
    "CoordinateSystem", {"spherical": 1, "cubed_sphere": 2, "cartesian": 3}
)
ZDirection = enum.Enum("ZDirection", {"sidechain": 1, "backbone": 2, "outward": 3})


def define_coordinate_system(pos_N, pos_CA, pos_C, z_direction):
    """Defines a local reference system based on N, CA, and C atom positions"""

    # Define local coordinate system
    e1 = pos_C - pos_N
    e1 /= np.linalg.norm(e1)

    # Define CB positions by rotating N atoms around CA-C axis 120 degr
    pos_N_res = pos_N - pos_CA
    axis = pos_CA - pos_C
    pos_CB = np.dot(
        Bio.PDB.rotaxis((120.0 / 180.0) * np.pi, Bio.PDB.vectors.Vector(axis)),
        pos_N_res,
    )
    e2 = pos_CB
    e2 /= np.linalg.norm(e2)
    e3 = np.cross(e1, e2)

    # N-C and e2 are not perfectly perpendical to one another. We adjust e2.
    e2 = np.cross(e1, -e3)

    if z_direction == ZDirection.outward:
        # Use e3 as z-direction
        rot_matrix = np.array([e1, e2, e3])
    elif z_direction == ZDirection.backbone:
        # Use backbone direction as z-direction
        rot_matrix = np.array([e2, e3, e1])
    elif z_direction == ZDirection.sidechain:
        # Use sidechain direction as z-direction
        rot_matrix = np.array([e3, e1, e2])
    else:
        raise "Unknown z-direction "

    return rot_matrix

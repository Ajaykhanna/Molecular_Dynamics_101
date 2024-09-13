#!/usr/bin/env python3

# ------------------------------------------------------------------------
# Python Script to Generate TeraChem Ground State Optimization
# Input Files in Explicit Solvent Environment
# For Multiple Dyes (>=2) in a sigle solvent
# By: Ajay Khanna | Sep.13.2024 [akhanna2@ucmerced.edu] | Isborn Lab
# ------------------------------------------------------------------------

import re
import os
import numpy as np
import argparse


# Function to parse command line arguments
def parse_args():
    """
    Parse command-line arguments for the script.

    Returns:
        argparse.Namespace: A namespace containing all the arguments.
    """
    parser = argparse.ArgumentParser(
        description="Extract Snapshots and Generate TeraChem Input Files for Multiple Dyes"
    )
    parser.add_argument(
        "--input",
        "-i",
        type=argparse.FileType("r"),
        required=True,
        help="MDTraj(xyz) file",
    )
    parser.add_argument(
        "--solv_charge",
        "-c",
        type=argparse.FileType("r"),
        required=True,
        help="Solvent point charge file (cols)",
    )
    parser.add_argument(
        "--qm_radius", "-r_qm", type=float, default=5, help="QM Radius (default=5.0A)"
    )
    parser.add_argument(
        "--mm_radius", "-r_mm", type=float, default=27, help="MM Radius (default=27.0A)"
    )
    parser.add_argument(
        "--nDyes", "-n_dyes", type=int, required=True, help="Number of Dyes"
    )
    parser.add_argument(
        "--dye_atoms",
        "-d_atoms",
        nargs="+",
        type=int,
        required=True,
        help="Total atoms in each dye (e.g., --dye_atoms 17 42)",
    )
    parser.add_argument(
        "--total_nDyes_atoms",
        "-tot_d_atoms",
        type=int,
        required=True,
        help="Total number of atoms in all dyes",
    )
    parser.add_argument(
        "--nAtoms_solvent",
        "-n_solvent_atoms",
        type=int,
        required=True,
        help="Number of atoms per solvent molecule",
    )
    parser.add_argument(
        "--total_frames",
        "-f",
        type=int,
        required=True,
        help="Total number of snapshots",
    )
    parser.add_argument(
        "--total_atoms", "-a", type=int, required=True, help="Total atoms in the box"
    )
    return parser.parse_args()


# Function to check if any atom in the solvent molecule is within qm_radius of any dye atom
def is_within_qm_radius(dye_coords_list, solvent_mol_coords, qm_radius):
    """
    Check if any atom in the solvent molecule is within the QM radius of any dye atom.

    Args:
        dye_coords_list (list of np.ndarray): List of coordinates for each dye's atoms.
        solvent_mol_coords (np.ndarray): Coordinates of the solvent molecule's atoms.
        qm_radius (float): QM radius in angstroms.

    Returns:
        bool: True if solvent molecule is within QM radius of any dye, False otherwise.
    """
    for dye_coords in dye_coords_list:
        distances = np.sqrt(
            np.sum(
                (dye_coords[:, np.newaxis, :] - solvent_mol_coords[np.newaxis, :, :])
                ** 2,
                axis=2,
            )
        )
        if np.any(distances < qm_radius):
            return True
    return False


# Function to prepend a line to a file
def prepend_line(filename, line):
    """
    Prepend a line to the beginning of a file.

    Args:
        filename (str): Path to the file.
        line (str): Line to prepend to the file.
    """
    with open(filename, "r+") as f:
        content = f.read()
        f.seek(0, 0)
        f.write(f"{line}\n\n{content}")


def remove_integers_from_symbol(atom_line: str) -> str:
    """
    The function `remove_integers_from_symbol` takes a string representing an atom line, removes any
    integers from the symbol part, and returns the modified atom line.

    Args:
        param atom_line: Atom line is a string that contains information about an atom
        :type atom_line: str
        :return: Returns a string where the first part has any integers removed
        and the rest of the parts are joined together with spaces.
    """
    parts = atom_line.split()
    symbol = re.sub(r"\d+", "", parts[0])
    return f"{symbol} {' '.join(parts[1:])}"


# Main logic to process snapshots
def process_snapshots(args):
    """
    Process each snapshot in the trajectory to generate QM and MM files, along with TeraChem input files.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    solvent_charge_list = np.loadtxt(args.solv_charge.name)
    total_atoms_snapshot = args.total_atoms
    lines = args.input.readlines()
    nDyes = args.nDyes
    dye_atoms_list = args.dye_atoms
    nAtoms_solvent = args.nAtoms_solvent

    if len(dye_atoms_list) != nDyes:
        raise ValueError(
            "Number of dyes does not match the number of dye_atoms provided."
        )

    for frame in range(args.total_frames):
        frame_dir = str(frame + 1)
        os.makedirs(frame_dir, exist_ok=True)

        qm_filename = os.path.join(
            frame_dir, f"solute_solvent_{args.qm_radius}ang_qm.xyz"
        )
        mm_filename = os.path.join(
            frame_dir, f"solute_solvent_{args.qm_radius}ang_mm.xyz"
        )
        tc_input_filename = os.path.join(
            frame_dir, f"tc_camb3lyp_{args.qm_radius}ang_opt_gs.in"
        )

        with open(tc_input_filename, "w") as params_file:
            params_file.write(terachem_params(mm_filename, qm_filename))

        dye_coords_list = []
        dye_atom_labels_list = []
        solvent_atom_labels = []
        solvent_coords = []

        # Calculate starting index for the current frame
        line_start = (
            frame * (total_atoms_snapshot + 2) + 2
        )  # Skip first two lines for each frame
        line_end = line_start + total_atoms_snapshot

        # Read dyes and solvents from the trajectory file
        idx = line_start

        # Read dye atoms
        for dye_idx in range(nDyes):
            dye_atom_count = dye_atoms_list[dye_idx]
            dye_coords = []
            dye_labels = []
            for _ in range(dye_atom_count):
                current_line = lines[idx].split()
                label = current_line[0]
                coords = [float(current_line[j]) for j in range(1, 4)]
                dye_coords.append(coords)
                dye_labels.append(label)
                idx += 1
            dye_coords_list.append(np.array(dye_coords))
            dye_atom_labels_list.append(dye_labels)

        # Read solvent atoms
        solvent_atom_count = total_atoms_snapshot - args.total_nDyes_atoms
        num_solvent_molecules = solvent_atom_count // nAtoms_solvent
        for _ in range(solvent_atom_count):
            current_line = lines[idx].split()
            label = current_line[0]
            coords = [float(current_line[j]) for j in range(1, 4)]
            solvent_coords.append(coords)
            solvent_atom_labels.append(label)
            idx += 1
        solvent_coords = np.array(solvent_coords)

        # Group solvent atoms into molecules
        solvent_molecules = []
        for i in range(0, len(solvent_coords), nAtoms_solvent):
            mol_coords = solvent_coords[i : i + nAtoms_solvent]
            mol_labels = solvent_atom_labels[i : i + nAtoms_solvent]
            solvent_molecules.append((mol_labels, mol_coords))

        # Determine which solvent molecules are within qm_radius of any dye
        qm_solvent_indices = []
        mm_solvent_indices = []
        for mol_idx, (mol_labels, mol_coords) in enumerate(solvent_molecules):
            if is_within_qm_radius(dye_coords_list, mol_coords, args.qm_radius):
                qm_solvent_indices.append(mol_idx)
            else:
                mm_solvent_indices.append(mol_idx)

        # Now write the QM and MM files, preserving the order
        with open(qm_filename, "w") as qm_file, open(mm_filename, "w") as mm_file:
            total_qm_atoms = 0

            # Write dye atoms to QM file
            for dye_labels, dye_coords in zip(dye_atom_labels_list, dye_coords_list):
                for label, coord in zip(dye_labels, dye_coords):
                    qm_file.write(
                        f"{remove_integers_from_symbol(label)}\t{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\n"
                    )
                    total_qm_atoms += 1

            # Write solvent molecules
            for mol_idx, (mol_labels, mol_coords) in enumerate(solvent_molecules):
                if mol_idx in qm_solvent_indices:
                    # Write solvent molecule to QM file
                    for label, coord in zip(mol_labels, mol_coords):
                        qm_file.write(
                            f"{remove_integers_from_symbol(label)}\t{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\n"
                        )
                        total_qm_atoms += 1
                else:
                    # Write solvent molecule to MM file with charges
                    for atom_idx_in_mol, coord in enumerate(mol_coords):
                        # Use modulo in case the charge list is shorter than the number of atoms per solvent molecule
                        charge = solvent_charge_list[
                            atom_idx_in_mol % len(solvent_charge_list)
                        ]
                        mm_file.write(
                            f"{charge:.6f}\t{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\n"
                        )

        # Prepend total number of atoms to QM file
        prepend_line(qm_filename, str(total_qm_atoms))

    print("------> Done!! <------")
    print("------> Buy Developer a Beer!! <------")


def terachem_params(point_charge_file, coordinate_file):
    """
    Generate TeraChem input parameters as a string.

    Args:
        point_charge_file (str): Filename of the MM point charges file.
        coordinate_file (str): Filename of the QM coordinates file.

    Returns:
        str: A string containing the TeraChem input parameters.
    """
    return f"""# TeraChem Job-Control Info.
charge        0
spinmult      1
basis         6-31G*
method        camb3lyp

# Type of Job Ground State Optimization
run           minimize
new_minimizer yes

# XYZ Filename/path
pointcharges  {os.path.basename(point_charge_file)}
coordinates   {os.path.basename(coordinate_file)}
dispersion    d3

# Scratch Directory Information
scrdir        gs_opt
keep_scr      yes

# Job Accuracy & Hardware Control
precision     double
threall       1.0e-15
convthre      3.0e-7
dftgrid       5
gpus          all
safemode      no
end
"""


if __name__ == "__main__":
    args = parse_args()
    process_snapshots(args)

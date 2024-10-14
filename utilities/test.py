#!/usr/bin/env python3

# ------------------------------------------------------------------------
# MD Trajectory to QM/MM Input Files Generator
# Updated to Utilize MM Radius
# By: [Your Name] | [Date] | [Your Email]
# ------------------------------------------------------------------------

import os
import numpy as np
import argparse


def parse_args():
    """
    Parses command-line arguments provided by the user.

    Returns:
        argparse.Namespace: An object containing all the parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Generate TeraChem and Gaussian Input Files for Multiple Dyes with Custom MM Radius"
    )
    parser.add_argument(
        "--input",
        "-i",
        type=argparse.FileType("r"),
        required=True,
        help="MDTraj (XYZ) file containing the trajectory",
    )
    parser.add_argument(
        "--solv_charge",
        "-c",
        type=argparse.FileType("r"),
        required=True,
        help="Solvent point charge file (one charge per line)",
    )
    parser.add_argument(
        "--qm_radius",
        "-r_qm",
        type=float,
        default=5,
        help="QM Radius in angstroms (default: 5.0 Å)",
    )
    parser.add_argument(
        "--mm_radius",
        "-r_mm",
        type=float,
        default=None,
        help="MM Radius in angstroms (default: include all solvent molecules outside QM region)",
    )
    parser.add_argument(
        "--nDyes",
        "-n_dyes",
        type=int,
        required=True,
        help="Number of dyes in the system",
    )
    parser.add_argument(
        "--dye_atoms",
        "-d_atoms",
        nargs="+",
        type=int,
        required=True,
        help="Number of atoms in each dye (e.g., --dye_atoms 17 42)",
    )
    parser.add_argument(
        "--total_nDyes_atoms",
        "-tot_d_atoms",
        type=int,
        required=True,
        help="Total number of atoms in all dyes combined",
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
        help="Total number of frames (snapshots) in the trajectory",
    )
    parser.add_argument(
        "--total_atoms",
        "-a",
        type=int,
        required=True,
        help="Total number of atoms per frame in the trajectory",
    )
    parser.add_argument(
        "--dye_MM_charges",
        nargs="*",
        help="List of dye_index:filename pairs to specify dyes to be converted to MM charges (e.g., 1:first_dye_MM_charge.txt)",
    )
    return parser.parse_args()


def is_within_radius(dye_coords_list, mol_coords, radius):
    """
    Checks if any atom in the molecule is within the given radius of any dye atom.

    Args:
        dye_coords_list (list of np.ndarray): List containing arrays of dye atom coordinates.
        mol_coords (np.ndarray): Coordinates of atoms in a molecule (solvent or dye).
        radius (float): The radius threshold.

    Returns:
        bool: True if the molecule is within the radius of any dye atom, False otherwise.
    """
    for dye_coords in dye_coords_list:
        # Compute distances between all dye atoms and molecule atoms
        distances = np.sqrt(
            np.sum(
                (dye_coords[:, np.newaxis, :] - mol_coords[np.newaxis, :, :]) ** 2,
                axis=2,
            )
        )
        if np.any(distances < radius):
            return True
    return False


def prepend_line(filename, line):
    """
    Prepends a line to the beginning of a file.

    Args:
        filename (str): Path to the file.
        line (str): Line to prepend.
    """
    with open(filename, "r+") as f:
        content = f.read()
        f.seek(0, 0)
        f.write(f"{line}\n\n{content}")


def process_snapshots(args):
    """
    Processes each frame in the trajectory, extracts dyes and solvent molecules within the QM and MM radii,
    and writes the QM and MM files while preserving the atom ordering. Also generates Gaussian input files.

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

    # Parse dye MM charges
    dye_MM_charge_files = (
        {}
    )  # Dictionary to store dye index and corresponding MM charge file
    if args.dye_MM_charges:
        for item in args.dye_MM_charges:
            try:
                index_str, filename = item.split(":")
                index = int(index_str)
                if index < 1 or index > nDyes:
                    raise ValueError(
                        f"Dye index {index} is out of range. Should be between 1 and {nDyes}."
                    )
                dye_MM_charge_files[index] = filename
            except ValueError:
                raise ValueError(
                    "Each item in --dye_MM_charges should be in the format dye_index:filename"
                )

    for frame in range(args.total_frames):
        frame_dir = str(frame + 1)
        os.makedirs(frame_dir, exist_ok=True)

        qm_filename = os.path.join(
            frame_dir, f"solute_solvent_{args.qm_radius}ang_qm_gs.xyz"
        )
        mm_filename = os.path.join(
            frame_dir, f"solute_solvent_{args.qm_radius}ang_mm.xyz"
        )
        tc_input_filename = os.path.join(
            frame_dir, f"tc_camb3lyp_{args.qm_radius}ang_opt_gs.in"
        )
        gaussian_input_filename = os.path.join(
            frame_dir, f"gaussian_input_{args.qm_radius}ang.com"
        )

        with open(tc_input_filename, "w") as params_file:
            params_file.write(terachem_params(mm_filename, qm_filename))

        dye_coords_list = []
        dye_atom_labels_list = []
        dye_atom_charges_list = []
        solvent_atom_labels = []
        solvent_coords = []

        # Calculate starting index for the current frame
        line_start = (
            frame * (total_atoms_snapshot + 2) + 2
        )  # Skip first two lines for each frame

        # Read dyes and solvents from the trajectory file
        idx = line_start

        # Read dye atoms
        for dye_idx in range(nDyes):
            dye_atom_count = dye_atoms_list[dye_idx]
            dye_coords = []
            dye_labels = []
            dye_charges = None
            dye_index = dye_idx + 1  # Dye indices start from 1

            if dye_index in dye_MM_charge_files:
                # Read charges from file
                dye_charge_file = dye_MM_charge_files[dye_index]
                dye_charges = np.loadtxt(dye_charge_file)
                if len(dye_charges) != dye_atom_count:
                    raise ValueError(
                        f"Number of charges in {dye_charge_file} does not match number of atoms in dye {dye_index}."
                    )

            for atom_idx in range(dye_atom_count):
                current_line = lines[idx].split()
                label = current_line[0]
                coords = [float(current_line[j]) for j in range(1, 4)]
                dye_coords.append(coords)
                if dye_charges is not None:
                    # Replace label with charge
                    charge = dye_charges[atom_idx]
                    label = f"{charge:.6f}"
                dye_labels.append(label)
                idx += 1

            dye_coords_list.append(np.array(dye_coords))
            dye_atom_labels_list.append(dye_labels)
            dye_atom_charges_list.append(dye_charges)

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

        # Determine which solvent molecules are within QM and MM radii
        # Exclude dyes converted to MM charges from QM and MM region determination
        qm_dye_coords_list = []
        mm_dye_coords_list = []
        for dye_index, (coords, charges) in enumerate(
            zip(dye_coords_list, dye_atom_charges_list), start=1
        ):
            if dye_index not in dye_MM_charge_files:
                qm_dye_coords_list.append(coords)
            else:
                mm_dye_coords_list.append(coords)

        # Solvent molecules within QM radius
        qm_solvent_indices = []
        # Solvent molecules within MM radius (but outside QM radius)
        mm_solvent_indices = []
        # Solvent molecules beyond MM radius
        discarded_solvent_indices = []

        for mol_idx, (mol_labels, mol_coords) in enumerate(solvent_molecules):
            if is_within_radius(qm_dye_coords_list, mol_coords, args.qm_radius):
                qm_solvent_indices.append(mol_idx)
            elif args.mm_radius is None or is_within_radius(
                qm_dye_coords_list, mol_coords, args.mm_radius
            ):
                mm_solvent_indices.append(mol_idx)
            else:
                discarded_solvent_indices.append(mol_idx)

        # Now write the QM and MM files, preserving the order
        with open(qm_filename, "w") as qm_file, open(mm_filename, "w") as mm_file:
            total_qm_atoms = 0

            # Write dye atoms
            for dye_index, (dye_labels, dye_coords, dye_charges) in enumerate(
                zip(dye_atom_labels_list, dye_coords_list, dye_atom_charges_list),
                start=1,
            ):
                if dye_index in dye_MM_charge_files:
                    # Check if dye is within MM radius
                    include_in_mm = True
                    if args.mm_radius is not None:
                        if not is_within_radius(
                            qm_dye_coords_list, dye_coords, args.mm_radius
                        ):
                            include_in_mm = False
                    if include_in_mm:
                        # Write to MM file with charges
                        for label, coord in zip(dye_labels, dye_coords):
                            mm_file.write(
                                f"{label}\t{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\n"
                            )
                else:
                    # Write to QM file
                    for label, coord in zip(dye_labels, dye_coords):
                        qm_file.write(
                            f"{label}\t{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\n"
                        )
                        total_qm_atoms += 1

            # Write solvent molecules
            for mol_idx, (mol_labels, mol_coords) in enumerate(solvent_molecules):
                if mol_idx in qm_solvent_indices:
                    # Write solvent molecule to QM file
                    for label, coord in zip(mol_labels, mol_coords):
                        qm_file.write(
                            f"{label}\t{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\n"
                        )
                        total_qm_atoms += 1
                elif mol_idx in mm_solvent_indices:
                    # Write solvent molecule to MM file with charges
                    for atom_idx_in_mol, coord in enumerate(mol_coords):
                        # Use modulo in case the charge list is shorter than the number of atoms per solvent molecule
                        charge = solvent_charge_list[
                            atom_idx_in_mol % len(solvent_charge_list)
                        ]
                        mm_file.write(
                            f"{charge:.6f}\t{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\n"
                        )
                # Else: discard solvent molecule (do not write it to any file)

        # Prepend total number of atoms to QM file
        prepend_line(qm_filename, str(total_qm_atoms))

        # Generate Gaussian input files
        generate_gaussian_input(
            gaussian_input_filename,
            dye_atom_labels_list,
            dye_coords_list,
            solvent_molecules,
            qm_solvent_indices,
            mm_solvent_indices,
            solvent_charge_list,
            dye_MM_charge_files,
        )

        # Generate Diabatization input files
        generate_diabatization_inputs(
            frame_dir,
            dye_atom_labels_list,
            dye_coords_list,
            solvent_molecules,
            mm_solvent_indices,
            solvent_charge_list,
        )

    print("------> Done!! <------")
    print("------> Buy Developer a Beer!! <------")


def terachem_params(point_charge_file, coordinate_file):
    """
    Generates the TeraChem input parameters as a formatted string.

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


def generate_gaussian_input(
    filename,
    dye_atom_labels_list,
    dye_coords_list,
    solvent_molecules,
    qm_solvent_indices,
    mm_solvent_indices,
    solvent_charge_list,
    dye_MM_charge_files,
):
    """
    Generates the Gaussian input file for a given frame.

    Args:
        filename (str): Path to the Gaussian input file to be created.
        dye_atom_labels_list (list): List of lists containing atom labels for each dye.
        dye_coords_list (list): List of numpy arrays containing coordinates for each dye.
        solvent_molecules (list): List of tuples containing solvent atom labels and coordinates.
        qm_solvent_indices (list): List of indices of solvent molecules in the QM region.
        mm_solvent_indices (list): List of indices of solvent molecules in the MM region.
        solvent_charge_list (np.ndarray): Array of solvent charge values.
        dye_MM_charge_files (dict): Dictionary mapping dye indices to MM charge files.
    """
    with open(filename, "w") as gauss_file:
        # Write checkpoint filename and job keywords
        chk_filename = os.path.splitext(os.path.basename(filename))[0] + ".chk"
        gauss_file.write(f"%chk={chk_filename}\n")
        gauss_file.write(
            "#p tda(nstates=6, root=1) cam-b3lyp/6-31g(d) nosymm charge\n\n"
        )
        gauss_file.write("Dyes in Solvent\n\n")
        gauss_file.write("0 1\n")

        # Write QM coordinates
        for dye_index, (dye_labels, dye_coords) in enumerate(
            zip(dye_atom_labels_list, dye_coords_list), start=1
        ):
            if dye_index not in dye_MM_charge_files:
                for label, coord in zip(dye_labels, dye_coords):
                    gauss_file.write(
                        f"{label}\t{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\n"
                    )

        for mol_idx in qm_solvent_indices:
            mol_labels, mol_coords = solvent_molecules[mol_idx]
            for label, coord in zip(mol_labels, mol_coords):
                gauss_file.write(
                    f"{label}\t{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\n"
                )

        gauss_file.write("\n")

        # Write MM coordinates with charges (X Y Z charge)
        # Include dyes converted to MM charges
        for dye_index, (dye_labels, dye_coords) in enumerate(
            zip(dye_atom_labels_list, dye_coords_list), start=1
        ):
            if dye_index in dye_MM_charge_files:
                # Check if dye is within MM radius
                include_in_mm = True
                if args.mm_radius is not None:
                    if not is_within_radius([dye_coords], dye_coords, args.mm_radius):
                        include_in_mm = False
                if include_in_mm:
                    for label, coord in zip(dye_labels, dye_coords):
                        # In label, we have the charge as a string
                        charge = float(label)
                        gauss_file.write(
                            f"{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\t{charge:.6f}\n"
                        )

        for mol_idx in mm_solvent_indices:
            mol_labels, mol_coords = solvent_molecules[mol_idx]
            for atom_idx_in_mol, coord in enumerate(mol_coords):
                # Use modulo in case the charge list is shorter than the number of atoms per solvent molecule
                charge = solvent_charge_list[atom_idx_in_mol % len(solvent_charge_list)]
                gauss_file.write(
                    f"{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\t{charge:.6f}\n"
                )

        gauss_file.write("\n")


def generate_diabatization_inputs(
    frame_dir,
    dye_atom_labels_list,
    dye_coords_list,
    solvent_molecules,
    mm_solvent_indices,
    solvent_charge_list,
):
    """
    Generates the Gaussian input files required for diabatization calculations.

    Args:
        frame_dir (str): Directory of the current frame.
        dye_atom_labels_list (list): List of lists containing atom labels for each dye.
        dye_coords_list (list): List of numpy arrays containing coordinates for each dye.
        solvent_molecules (list): List of solvent molecules (labels and coordinates).
        mm_solvent_indices (list): List of indices of solvent molecules in the MM region.
        solvent_charge_list (np.ndarray): Array of solvent charge values.
    """
    # Common Gaussian keywords
    gaussian_keywords = "#p TDA(nstates=6, root=1) charge cam-b3lyp/6-31g(d) nosymm EmpiricalDispersion=GD3"

    # File 1: All dyes with QM and MM coordinates
    filename_all = os.path.join(frame_dir, "diabat_all_dyes.com")
    with open(filename_all, "w") as gauss_file:
        write_gaussian_header(gauss_file, filename_all, gaussian_keywords)
        gauss_file.write("0 1\n")

        # Write all dyes
        for dye_labels, dye_coords in zip(dye_atom_labels_list, dye_coords_list):
            for label, coord in zip(dye_labels, dye_coords):
                gauss_file.write(
                    f"{label}\t{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\n"
                )

        gauss_file.write("\n")

        # Write MM coordinates
        write_mm_coordinates(
            gauss_file,
            [],
            [],
            solvent_molecules,
            mm_solvent_indices,
            solvent_charge_list,
            zero_charges=False,
        )

    # File 2: Only first dye, other dyes as point charges with zero charges
    filename_dye1 = os.path.join(frame_dir, "diabat_dye1.com")
    with open(filename_dye1, "w") as gauss_file:
        write_gaussian_header(gauss_file, filename_dye1, gaussian_keywords)
        gauss_file.write("0 1\n")

        # Write first dye
        dye_labels = dye_atom_labels_list[0]
        dye_coords = dye_coords_list[0]
        for label, coord in zip(dye_labels, dye_coords):
            gauss_file.write(
                f"{label}\t{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\n"
            )

        gauss_file.write("\n")

        # Write MM coordinates
        # Include other dyes as point charges with zero charges
        other_dyes_labels = dye_atom_labels_list[1:]
        other_dyes_coords = dye_coords_list[1:]
        write_mm_coordinates(
            gauss_file,
            other_dyes_labels,
            other_dyes_coords,
            solvent_molecules,
            mm_solvent_indices,
            solvent_charge_list,
            zero_charges=True,
        )

    # File 3: Only second dye, other dyes as point charges with zero charges
    if len(dye_atom_labels_list) > 1:
        filename_dye2 = os.path.join(frame_dir, "diabat_dye2.com")
        with open(filename_dye2, "w") as gauss_file:
            write_gaussian_header(gauss_file, filename_dye2, gaussian_keywords)
            gauss_file.write("0 1\n")

            # Write second dye
            dye_labels = dye_atom_labels_list[1]
            dye_coords = dye_coords_list[1]
            for label, coord in zip(dye_labels, dye_coords):
                gauss_file.write(
                    f"{label}\t{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\n"
                )

            gauss_file.write("\n")

            # Write MM coordinates
            # Include other dyes as point charges with zero charges
            other_dyes_labels = dye_atom_labels_list[:1] + dye_atom_labels_list[2:]
            other_dyes_coords = dye_coords_list[:1] + dye_coords_list[2:]
            write_mm_coordinates(
                gauss_file,
                other_dyes_labels,
                other_dyes_coords,
                solvent_molecules,
                mm_solvent_indices,
                solvent_charge_list,
                zero_charges=True,
            )


def write_gaussian_header(gauss_file, filename, keywords):
    """
    Writes the header for Gaussian input files.

    Args:
        gauss_file (file object): Open file object for writing.
        filename (str): Name of the Gaussian input file.
        keywords (str): Gaussian job keywords.
    """
    chk_filename = os.path.splitext(os.path.basename(filename))[0] + ".chk"
    gauss_file.write(f"%chk={chk_filename}\n")
    gauss_file.write(f"{keywords}\n\n")
    gauss_file.write("Dyes in Solvent\n\n")


def write_mm_coordinates(
    gauss_file,
    dye_atom_labels_list,
    dye_coords_list,
    solvent_molecules,
    mm_solvent_indices,
    solvent_charge_list,
    zero_charges=False,
):
    """
    Writes MM coordinates with charges into the Gaussian input file.

    Args:
        gauss_file (file object): Open file object for writing.
        dye_atom_labels_list (list): List of lists of dye atom labels.
        dye_coords_list (list): List of numpy arrays of dye coordinates.
        solvent_molecules (list): List of solvent molecules (labels and coordinates).
        mm_solvent_indices (list): List of indices of solvent molecules in the MM region.
        solvent_charge_list (np.ndarray): Array of solvent charge values.
        zero_charges (bool): If True, set charges to zero for dyes.
    """
    # Write dyes as MM charges
    for dye_labels, dye_coords in zip(dye_atom_labels_list, dye_coords_list):
        for coord in dye_coords:
            charge = 0.0 if zero_charges else None
            if charge is not None:
                gauss_file.write(
                    f"{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\t{charge:.6f}\n"
                )

    # Write solvent molecules
    for mol_idx in mm_solvent_indices:
        mol_labels, mol_coords = solvent_molecules[mol_idx]
        for atom_idx_in_mol, coord in enumerate(mol_coords):
            # Use modulo in case the charge list is shorter than the number of atoms per solvent molecule
            charge = solvent_charge_list[atom_idx_in_mol % len(solvent_charge_list)]
            gauss_file.write(
                f"{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\t{charge:.6f}\n"
            )
    gauss_file.write("\n")


if __name__ == "__main__":
    args = parse_args()
    process_snapshots(args)

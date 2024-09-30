import re
import os
import numpy as np


def dummy_gaussian_input_file(filename: str):
    """
    Generates a dummy Gaussian input file for testing purposes.
    """

    with open(filename, "w") as f:
        f.write(
            """
%chk=dummy.chk
#p HF/6-31g* nosymm
                
Dummy Gaussian Input File
                
0 1
C 0.000000 0.000000 0.000000
H 0.000000 0.000000 1.089000
H 1.026719 0.000000 -0.363000
H -0.513360 -0.889165 -0.363000
H -0.513360 0.889165 -0.36300
                

"""
        )


def remove_integers_from_symbol(atom_line: str) -> str:
    """
    Removes any integers from the symbol part of an atom line.

    Args:
        atom_line (str): String representing an atom line.
    Returns:
        str: Atom symbol without integers.
    """
    symbol = re.sub(r"\d+", "", atom_line.split()[0])
    return symbol


def write_qm_coordinates(
    gauss_file,
    dye_atom_labels_list: list[list[str]],
    dye_coords_list: list[np.ndarray],
    solvent_molecules: list[tuple[list[str], np.ndarray]],
    qm_solvent_indices: list[int],
):
    """
    Writes the QM coordinates to the Gaussian input file.

    Args:
        gauss_file: File object to write to.
        dye_atom_labels_list: List of lists containing atom labels for each dye.
        dye_coords_list: List of numpy arrays containing coordinates for each dye.
        solvent_molecules: List of tuples containing solvent atom labels and coordinates.
        qm_solvent_indices: List of indices of solvent molecules in the QM region.
    """
    for dye_labels, dye_coords in zip(dye_atom_labels_list, dye_coords_list):
        for label, coord in zip(dye_labels, dye_coords):
            gauss_file.write(
                f"{remove_integers_from_symbol(label)}\t"
                f"{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\n"
            )

    for mol_idx in qm_solvent_indices:
        mol_labels, mol_coords = solvent_molecules[mol_idx]
        for label, coord in zip(mol_labels, mol_coords):
            gauss_file.write(
                f"{remove_integers_from_symbol(label)}\t"
                f"{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\n"
            )


def write_mm_coordinates(
    gauss_file,
    solvent_molecules: list[tuple[list[str], np.ndarray]],
    mm_solvent_indices: list[int],
    solvent_charge_list: np.ndarray,
):
    """
    Writes the MM coordinates with charges to the Gaussian input file.

    Args:
        gauss_file: File object to write to.
        solvent_molecules: List of tuples containing solvent atom labels and coordinates.
        mm_solvent_indices: List of indices of solvent molecules in the MM region.
        solvent_charge_list: Array of solvent charge values.
    """
    for mol_idx in mm_solvent_indices:
        mol_labels, mol_coords = solvent_molecules[mol_idx]
        for atom_idx_in_mol, coord in enumerate(mol_coords):
            charge = solvent_charge_list[atom_idx_in_mol % len(solvent_charge_list)]
            gauss_file.write(
                f"{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\t{charge:.6f}\n"
            )


def generate_gaussian_input_file(
    filename: str,
    dye_atom_labels_list: list[list[str]],
    dye_coords_list: list[np.ndarray],
    solvent_molecules: list[tuple[list[str], np.ndarray]],
    qm_solvent_indices: list[int],
    mm_solvent_indices: list[int],
    solvent_charge_list: np.ndarray,
    net_charge: int = 0,
    spin_mult: int = 1,
    route_section: str = "",
    title: str = "",
    header_options: str = "",
    extra_sections: str = "",
):
    """
    Generic function to generate Gaussian input files.

    Args:
        filename (str): Path to the Gaussian input file to be created.
        dye_atom_labels_list (list): List of lists containing atom labels for each dye.
        dye_coords_list (list): List of numpy arrays containing coordinates for each dye.
        solvent_molecules (list): List of tuples containing solvent atom labels and coordinates.
        qm_solvent_indices (list): List of indices of solvent molecules in the QM region.
        mm_solvent_indices (list): List of indices of solvent molecules in the MM region.
        solvent_charge_list (np.ndarray): Array of solvent charge values.
        net_charge (int): Net charge of the system.
        spin_mult (int): Spin multiplicity of the system.
        route_section (str): Gaussian route section (without the #p).
        title (str): Title section of the Gaussian input file.
        header_options (str): Additional header options (e.g., %chk=filename).
        extra_sections (str): Any extra sections to add at the end of the input file.
    """
    with open(filename, "w") as gauss_file:
        # Write header
        gauss_file.write(f"{header_options}\n")
        gauss_file.write(f"#p {route_section}\n\n")
        gauss_file.write(f"{title}\n\n")
        gauss_file.write(f"{net_charge} {spin_mult}\n")

        # Write QM coordinates
        write_qm_coordinates(
            gauss_file,
            dye_atom_labels_list,
            dye_coords_list,
            solvent_molecules,
            qm_solvent_indices,
        )

        gauss_file.write("\n")

        # Write MM coordinates with charges
        write_mm_coordinates(
            gauss_file,
            solvent_molecules,
            mm_solvent_indices,
            solvent_charge_list,
        )

        gauss_file.write("\n")

        if extra_sections:
            gauss_file.write(extra_sections)
            gauss_file.write("\n")


def generate_ground_state_energy_files(
    filename: str,
    dye_atom_labels_list: list[list[str]],
    dye_coords_list: list[np.ndarray],
    solvent_molecules: list[tuple[list[str], np.ndarray]],
    qm_solvent_indices: list[int],
    mm_solvent_indices: list[int],
    solvent_charge_list: np.ndarray,
    net_charge: int = 0,
    spin_mult: int = 1,
    dft_func: str = "cam-b3lyp",
    basis: str = "6-31g*",
):
    """
    Generates the Gaussian input file for ground state energy calculations.
    """
    chk_filename = os.path.splitext(os.path.basename(filename))[0] + ".chk"
    header_options = f"%chk={chk_filename}"
    route_section = f"{dft_func}/{basis} nosymm charge EmpiricalDispersion=GD3"
    title = "Ground State Energy Calculations"

    generate_gaussian_input_file(
        filename,
        dye_atom_labels_list,
        dye_coords_list,
        solvent_molecules,
        qm_solvent_indices,
        mm_solvent_indices,
        solvent_charge_list,
        net_charge,
        spin_mult,
        route_section,
        title,
        header_options,
    )


def generate_ground_state_optimization_files(
    filename: str,
    dye_atom_labels_list: list[list[str]],
    dye_coords_list: list[np.ndarray],
    solvent_molecules: list[tuple[list[str], np.ndarray]],
    qm_solvent_indices: list[int],
    mm_solvent_indices: list[int],
    solvent_charge_list: np.ndarray,
    net_charge: int = 0,
    spin_mult: int = 1,
    dft_func: str = "cam-b3lyp",
    basis: str = "6-31g*",
    opt_freq: bool = False,
):
    """
    Generates the Gaussian input file for ground state optimization calculations.
    """
    chk_filename = os.path.splitext(os.path.basename(filename))[0] + ".chk"
    header_options = f"%chk={chk_filename}"
    if opt_freq:
        route_section = f"opt freq {dft_func}/{basis} nosymm charge EmpiricalDispersion=GD3"
        title = "Ground State Optimization and Frequency Calculations"
    else:
        route_section = f"opt {dft_func}/{basis} nosymm charge EmpiricalDispersion=GD3"
        title = "Ground State Optimization Calculations"

    generate_gaussian_input_file(
        filename,
        dye_atom_labels_list,
        dye_coords_list,
        solvent_molecules,
        qm_solvent_indices,
        mm_solvent_indices,
        solvent_charge_list,
        net_charge,
        spin_mult,
        route_section,
        title,
        header_options,
    )


def generate_ground_state_frequency_files(
    filename: str,
    dye_atom_labels_list: list[list[str]],
    dye_coords_list: list[np.ndarray],
    solvent_molecules: list[tuple[list[str], np.ndarray]],
    qm_solvent_indices: list[int],
    mm_solvent_indices: list[int],
    solvent_charge_list: np.ndarray,
    net_charge: int = 0,
    spin_mult: int = 1,
    dft_func: str = "cam-b3lyp",
    basis: str = "6-31g*",
):
    """
    Generates the Gaussian input file for ground state frequency calculations.
    """
    chk_filename = os.path.splitext(os.path.basename(filename))[0] + ".chk"
    header_options = f"%chk={chk_filename}"
    route_section = f"freq=(saveNM, HPModes) {dft_func}/{basis} nosymm charge EmpiricalDispersion=GD3"
    title = "Ground State Frequency Calculations"

    generate_gaussian_input_file(
        filename,
        dye_atom_labels_list,
        dye_coords_list,
        solvent_molecules,
        qm_solvent_indices,
        mm_solvent_indices,
        solvent_charge_list,
        net_charge,
        spin_mult,
        route_section,
        title,
        header_options,
    )


def generate_vertical_excitation_energy_file(
    filename: str,
    dye_atom_labels_list: list[list[str]],
    dye_coords_list: list[np.ndarray],
    solvent_molecules: list[tuple[list[str], np.ndarray]],
    qm_solvent_indices: list[int],
    mm_solvent_indices: list[int],
    solvent_charge_list: np.ndarray,
    net_charge: int = 0,
    spin_mult: int = 1,
    dft_func: str = "cam-b3lyp",
    basis: str = "6-31g*",
    nstates: int = 6,
    root: int = 1,
):
    """
    Generates the Gaussian input file for vertical excitation energy calculations.
    """
    chk_filename = os.path.splitext(os.path.basename(filename))[0] + ".chk"
    header_options = f"%chk={chk_filename}"
    route_section = (
        f"tda(nstates={nstates}, root={root}) {dft_func}/{basis} nosymm charge EmpiricalDispersion=GD3"
    )
    title = "Vertical Excitaion Energy Calculations"

    generate_gaussian_input_file(
        filename,
        dye_atom_labels_list,
        dye_coords_list,
        solvent_molecules,
        qm_solvent_indices,
        mm_solvent_indices,
        solvent_charge_list,
        net_charge,
        spin_mult,
        route_section,
        title,
        header_options,
    )


def generate_excited_state_optimization_files(
    filename: str,
    dye_atom_labels_list: list[list[str]],
    dye_coords_list: list[np.ndarray],
    solvent_molecules: list[tuple[list[str], np.ndarray]],
    qm_solvent_indices: list[int],
    mm_solvent_indices: list[int],
    solvent_charge_list: np.ndarray,
    net_charge: int = 0,
    spin_mult: int = 1,
    dft_func: str = "cam-b3lyp",
    basis: str = "6-31g*",
    nstates: int = 6,
    root: int = 1,
    opt_freq: bool = False,
):
    """
    Generates the Gaussian input file for excited state optimization calculations.
    """
    chk_filename = os.path.splitext(os.path.basename(filename))[0] + ".chk"
    header_options = f"%chk={chk_filename}"
    if opt_freq:
        route_section = f"tda(nstates={nstates}, root={root}) opt freq(saveNM, HPModes) {dft_func}/{basis} nosymm charge EmpiricalDispersion=GD3"
        title = "Excited State Optimization and Frequency Calculations"
    else:
        route_section = (
            f"tda(nstates={nstates}, root={root}) opt {dft_func}/{basis} nosymm charge EmpiricalDispersion=GD3"
        )
        title = "Excited State Optimization Calculations"

    generate_gaussian_input_file(
        filename,
        dye_atom_labels_list,
        dye_coords_list,
        solvent_molecules,
        qm_solvent_indices,
        mm_solvent_indices,
        solvent_charge_list,
        net_charge,
        spin_mult,
        route_section,
        title,
        header_options,
    )


def generate_excited_state_frequency_files(
    filename: str,
    dye_atom_labels_list: list[list[str]],
    dye_coords_list: list[np.ndarray],
    solvent_molecules: list[tuple[list[str], np.ndarray]],
    qm_solvent_indices: list[int],
    mm_solvent_indices: list[int],
    solvent_charge_list: np.ndarray,
    net_charge: int = 0,
    spin_mult: int = 1,
    dft_func: str = "cam-b3lyp",
    basis: str = "6-31g*",
    nstates: int = 6,
    root: int = 1,
):
    """
    Generates the Gaussian input file for excited state frequency calculations.
    """
    chk_filename = os.path.splitext(os.path.basename(filename))[0] + ".chk"
    header_options = f"%chk={chk_filename}"
    route_section = (
        f"freq=(saveNM, HPModes) tda(nstates={nstates}, root={root})\n"
        f"{dft_func}/{basis} nosymm charge EmpiricalDispersion=GD3"
    )
    title = "Excited State Frequency Calculations"

    generate_gaussian_input_file(
        filename,
        dye_atom_labels_list,
        dye_coords_list,
        solvent_molecules,
        qm_solvent_indices,
        mm_solvent_indices,
        solvent_charge_list,
        net_charge,
        spin_mult,
        route_section,
        title,
        header_options,
    )


def generate_charge_files(
    filename: str,
    dye_atom_labels_list: list[list[str]],
    dye_coords_list: list[np.ndarray],
    solvent_molecules: list[tuple[list[str], np.ndarray]],
    qm_solvent_indices: list[int],
    mm_solvent_indices: list[int],
    solvent_charge_list: np.ndarray,
    net_charge: int = 0,
    spin_mult: int = 1,
    theory: str = "HF",
    basis: str = "6-31g*",
    method: str = "MK",
    excited_state: bool = False,
    nstates: int = 6,
    root: int = 1,
):
    """
    Generates the Gaussian input file for charge calculations.
    """
    chk_filename = os.path.splitext(os.path.basename(filename))[0] + ".chk"
    header_options = f"%chk={chk_filename}"

    if method == "MK":
        pop_keyword = "(regular, MK)"
    elif method == "ChelpG":
        pop_keyword = "(regular, ChelpG)"
    else:
        pop_keyword = "(SaveNTO)"

    scf_options = "SCF=tight test"
    iop_options = "iop(6/33=2) iop(6/42=6) iop(6/50=1)"

    if excited_state:
        route_section = (
            f"tda(nstates={nstates}, root={root}) pop={pop_keyword}\n"
            f"{theory}/{basis} nosymm charge {iop_options} {scf_options} EmpiricalDispersion=GD3"
        )
        title = "Excited State Charge Calculations"
    else:
        route_section = (
            f"pop={pop_keyword} {theory}/{basis} nosymm charge "
            f"{iop_options} {scf_options}"
        )
        title = "Ground State Charge Calculations"

    generate_gaussian_input_file(
        filename,
        dye_atom_labels_list,
        dye_coords_list,
        solvent_molecules,
        qm_solvent_indices,
        mm_solvent_indices,
        solvent_charge_list,
        net_charge,
        spin_mult,
        route_section,
        title,
        header_options,
    )


def generate_transition_charge_files(
    filename: str,
    dye_atom_labels_list: list[list[str]],
    dye_coords_list: list[np.ndarray],
    solvent_molecules: list[tuple[list[str], np.ndarray]],
    qm_solvent_indices: list[int],
    mm_solvent_indices: list[int],
    solvent_charge_list: np.ndarray,
    net_charge: int = 0,
    spin_mult: int = 1,
    theory: str = "HF",
    basis: str = "6-31g*",
    method: str = "SaveNTO",
    nstates: int = 6,
    root: int = 1,
):
    """
    Generates the Gaussian input file for charge calculations.
    """
    chk_filename = os.path.splitext(os.path.basename(filename))[0] + ".chk"
    header_options = f"%chk={chk_filename}"

    if method == "MK":
        pop_keyword = "(regular, MK)"
    elif method == "ChelpG":
        pop_keyword = "(regular, ChelpG)"
    else:
        pop_keyword = f"SaveNTO density(transition={root})"

    route_section = (
        f"tda(nstates={nstates}, root={root}) pop={pop_keyword}\n"
        f"{theory}/{basis} nosymm charge EmpiricalDispersion=GD3"
    )
    title = "Transition State Charge Calculations"

    generate_gaussian_input_file(
        filename,
        dye_atom_labels_list,
        dye_coords_list,
        solvent_molecules,
        qm_solvent_indices,
        mm_solvent_indices,
        solvent_charge_list,
        net_charge,
        spin_mult,
        route_section,
        title,
        header_options,
    )


if __name__ == "__main__":
    print("Script to Generate Gaussian Input Files")
    dummy_gaussian_input_file("dummy.com")

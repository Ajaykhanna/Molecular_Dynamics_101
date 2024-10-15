# MD Trajectory to QM/MM Input Files Generator

## Overview

**MD Trajectory to QM/MM Input Files Generator** is a Python script designed to process molecular dynamics (MD) trajectory data containing multiple dyes and solvent molecules. It automates the extraction of dyes and nearby solvent molecules within a specified radius and generates input files for TeraChem and Gaussian software, including files necessary for diabatization calculations. The script maintains the original atom ordering and offers flexibility in configuring dyes as MM point charges.

---

## Features

- **Supports Multiple Dyes**: Handles systems with any number of dyes, each with a specified number of atoms.
- **QM and MM Region Determination**: Automatically assigns atoms to QM or MM regions based on distance criteria.
- **Flexible Dye Configuration**: Allows users to convert any dye into MM point charges by providing charge files.
- **Generates Input Files**:
  - **TeraChem Input Files**: QM and MM coordinates, along with TeraChem job control parameters.
  - **Gaussian Input Files**: Includes both QM and MM coordinates, formatted correctly for Gaussian simulations.
  - **Diabatization Input Files**: Prepares specialized Gaussian input files for diabatization calculations.
- **Maintains Atom Ordering**: Preserves the original ordering of atoms from the trajectory in all output files.
- **User-Friendly Command-Line Interface**: Provides clear and customizable command-line arguments.

---

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
  - [Command-Line Arguments](#command-line-arguments)
  - [Example](#example)
- [Input File Formats](#input-file-formats)
  - [Trajectory File](#trajectory-file)
  - [Solvent Charge File](#solvent-charge-file)
  - [Dye MM Charge Files](#dye-mm-charge-files)
- [Output Files](#output-files)
  - [TeraChem Files](#terachem-files)
  - [Gaussian Files](#gaussian-files)
  - [Diabatization Files](#diabatization-files)
- [Notes](#notes)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

---

## Installation

1. **Clone the Repository**:

   ```bash
   git clone https://github.com/yourusername/md-trajectory-to-qmmm-input.git
   ```

2. **Navigate to the Directory**:

   ```bash
   cd md-trajectory-to-qmmm-input
   ```

3. **Install Dependencies**:

   The script requires Python 3 and the following Python libraries:

   - `numpy`

   Install the dependencies using `pip`:

   ```bash
   pip install numpy
   ```

---

## Usage

The script is executed via the command line and requires several arguments to function correctly.

### Command-Line Arguments

- `--input` (`-i`): **(Required)** The trajectory file in XYZ format containing the frames to process. This XYZ can have either alphanumeric atomic symbols (C1 X Y Z) or just atomic sybmols (C X Y Z). For Alphanumeric XYZ format, the script will automatically convert to atomic symbols when passed "--xyz_alphanumeric or --xyzint" argument.
- `--solv_charge` (`-c`): **(Required)** File containing solvent point charge values (one per line).
- `--qm_radius` (`-r_qm`): **(Optional)** The radius in angstroms to include solvent molecules in the QM region around each dye (default: `5.0` Å).
- `--mm_radius` (`-r_mm`): **(Optional)** The radius for the MM region (default: None (Include all solvent molecules)).
- `--nDyes` (`-n_dyes`): **(Required)** Number of dyes present in the trajectory.
- `--dye_atoms` (`-d_atoms`): **(Required)** Number of atoms in each dye. Provide one integer per dye (e.g., `--dye_atoms 17 42`).
- `--total_nDyes_atoms` (`-tot_d_atoms`): **(Required)** Total number of atoms in all dyes combined.
- `--nAtoms_solvent` (`-n_solvent_atoms`): **(Required)** Number of atoms per solvent molecule.
- `--total_frames` (`-f`): **(Required)** Total number of frames (snapshots) in the trajectory file to process.
- `--total_atoms` (`-a`): **(Required)** Total number of atoms per frame in the trajectory file.
- `--dye_MM_charges`: **(Optional)** List of `dye_index:filename` pairs to specify dyes to be converted to MM charges (e.g., `1:first_dye_MM_charge.txt`).
- `--net_charge` (`-q`): **(Optional)** Net charge of the system (default: `0`).
- `--spin_mult` (`-s`): **(Optional)** Spin multiplicity of the system (default: `1`).
- `--gaussian_inputs` (-gau_inputs): **(Optional)** Generate Gaussian input files for various QM calculations (default = False).
- `--teracheem_inputs` (-tera_inputs): **(Optional)** Generate TeraChem input files for various QM calculations (default = False).
- `--xyz_alphanumeric` (-xyzint): **(Optional)** Use alphanumeric atomic symbols (C1 X Y Z) instead of atomic symbols (C X Y Z) in the XYZ file (default = False).

### Example

#### With one dye/ligand/molecule/choromophore

```bash
python md_to_qmmm_input.py --input twoframes_example.xyz --solv_charge solv_charge.txt \
--qm_radius 5 --nDyes 1 --dye_atoms 17 --total_nDyes_atoms 17 \
--nAtoms_solvent 10 --total_frames 10 --total_atoms 4489
```

#### With two or more dye/ligand/molecule/choromophore

```bash
python md_to_qmmm_input.py --input twoframes_example.xyz --solv_charge solv_charge.txt \
--qm_radius 5 --nDyes 2 --dye_atoms 17 42 --total_nDyes_atoms 59 \
--nAtoms_solvent 10 --total_frames 10 --total_atoms 4489
```

#### Defining a Custom MM Radius

```bash
python md_to_qmmm_input.py --input twoframes_example.xyz --solv_charge solv_charge.txt \
--qm_radius 5 --mm_radius 10 --nDyes 2 --dye_atoms 17 42 \
--total_nDyes_atoms 59 --nAtoms_solvent 3 --total_frames 10 \
--total_atoms 1000
```

#### Converting a Dye to MM Charges

```bash
python md_to_qmmm_input.py --input twoframes_example.xyz --solv_charge solv_charge.txt \
--qm_radius 5 --nDyes 2 --dye_atoms 17 42 --total_nDyes_atoms 59 \
--nAtoms_solvent 3 --total_frames 10 --total_atoms 4489 \
--dye_MM_charges 1:first_dye_MM_charge.txt --net_charge 0 --spin_mult 1
```

#### Generating Diabatization Input Files

```bash
python md_to_qmmm_input.py --input twoframes_example.xyz --solv_charge solv_charge.txt \
--qm_radius 5 --nDyes 2 --dye_atoms 17 42 --total_nDyes_atoms 59 \
--nAtoms_solvent 3 --total_frames 10 --total_atoms 4489 --net_charge 0 --spin_mult 1\
--gaussian_inputs
```

## Note

Currently, ``--dye_MM_charges`` and ``--gaussian_inputs`` or ``--terachem_inputs`` are not compatible. Please use one at a time.

---

## Input File Formats

### Trajectory File

- **Format**: XYZ format with a specific structure per frame.
- **Per Frame**:
  - **First Line**: Total number of atoms per frame (integer).
  - **Second Line**: Comment or empty line.
  - **Subsequent Lines**: Atom label and XYZ coordinates, with dyes first and then solvents.

**Example**:

```
1000
Frame 1
C1    0.000000    0.000000    0.000000
H1    0.629118    0.629118    0.629118
...
O1    1.234567    2.345678    3.456789
```

### Solvent Charge File

- **Content**: Charge values corresponding to each atom in a single solvent molecule.
- **Format**: One charge value per line.

**Example** (for a solvent molecule with 3 atoms):

```
0.423
-0.846
0.423
```

### Dye MM Charge Files

- **Purpose**: Contains charge values for each atom in a dye when converting the dye to MM point charges.
- **Format**: One charge value per line, matching the number of atoms in the dye.

---

## Output Files

For each frame processed, the script generates several files organized within a directory named after the frame number (`1`, `2`, ..., `N`).

### TeraChem Files

1. **QM XYZ File** (`solute_solvent_<qm_radius>ang_qm_gs.xyz`):

   - Contains the QM region atoms: dyes (unless converted to MM charges) and solvent molecules within the QM radius.
   - Atom labels and coordinates are preserved from the input trajectory.

2. **MM Point Charges File** (`solute_solvent_<qm_radius>ang_mm.xyz`):

   - Contains the MM region atoms: dyes converted to MM charges and solvent molecules outside the QM radius.
   - Atom labels are replaced with charge values.

3. **TeraChem Input File** (`tc_camb3lyp_<qm_radius>ang_opt_gs.in`):

   - Contains job control parameters for TeraChem.
   - References the QM and MM files generated.

### Gaussian Files

1. **Standard Gaussian Input File** (`gaussian_input_<qm_radius>ang.com`):

   - Includes both QM and MM coordinates.
   - MM coordinates are formatted as `X Y Z charge`.
   - Contains a header with the checkpoint filename and job keywords.

### Diabatization Files

1. **All Dyes File** (`diabat_all_dyes.com`):

   - Contains all dyes with QM and MM coordinates.
   - Uses specified Gaussian keywords for diabatization.

2. **First Dye File** (`diabat_dye1.com`):

   - Contains only the first dye in the QM region.
   - Other dyes are included as point charges with zero charges.

3. **Second Dye File** (`diabat_dye2.com`):

   - Contains only the second dye in the QM region.
   - Other dyes are included as point charges with zero charges.

---

## Notes

- **Atom Ordering**: The script maintains the ordering of atoms as they appear in the input trajectory file.
- **Dye Indices**: Dyes are indexed starting from 1.
- **Charge Files**: When converting dyes to MM charges, ensure the charge files have one charge per atom in the dye.
- **Solvent Charge List**: The solvent charge file should contain charges for one solvent molecule, repeated as necessary.
- **Trajectory Consistency**: Ensure that the `--total_atoms` argument matches the number of atoms specified in the trajectory file's first line per frame.

---

## Contributing

Contributions are welcome! Please follow these steps:

1. **Fork the Repository**: Click the "Fork" button at the top right of this page.
2. **Clone Your Fork**:

   ```bash
   git clone https://github.com/Ajaykhanna/Molecular_Dynamics_101.git
   ```

3. **Create a Branch**:

   ```bash
   git checkout -b feature/review-changes
   ```

4. **Make Changes**: Implement your feature or fix.
5. **Commit Changes**:

   ```bash
   git commit -am 'Add new feature'
   ```

6. **Push to Your Fork**:

   ```bash
   git push origin feature/review-changes
   ```

7. **Submit a Pull Request**: Go to the original repository and create a pull request.

---

## License

This project is licensed under the [MIT License](LICENSE).

---

## Contact

- **Author**: [Ajay Khanna]
- **Email**: [akhanna2@ucmerced.edu]
- **GitHub**: [Ajaykhanna](https://github.com/Ajaykhanna)

For any questions or issues, please open an issue on GitHub or contact me via email.

---

**Acknowledgments**:

- Thank you to all contributors and users who have provided feedback and suggestions.
- This script was inspired by the need to automate the preparation of QM/MM simulations in computational chemistry.

---

**Disclaimer**:

- Ensure that all input files are correctly formatted and that the command-line arguments match your system's specifications.
- The script is provided "as-is" without warranty of any kind.

---
**Happy QM/MM!**

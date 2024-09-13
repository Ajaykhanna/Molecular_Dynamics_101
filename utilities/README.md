## Extracting Snapshots in QM and MM region as XYZ files

___

This script processes a molecular dynamics trajectory file containing multiple dyes and solvent molecules. It extracts the dyes and solvent molecules within a specified QM radius and generates QM and MM files for TeraChem calculations.

How to Run the Script:

```python
python extract_qm_mm_snapshots.py --input trajectory.xyz --solv_charge charges.txt --qm_radius 5.0 --mm_radius 27.0 --nDyes 2 --dye_atoms 17 42 --total_nDyes_atoms 59 --nAtoms_solvent 10 --total_frames 2 --total_atoms 4489
```

### Prerequisites

* Python 3.x
* NumPy library installed
* An XYZ trajectory file in the specified format
* A solvent charge file containing charge values

### Arguments

```python
--input (-i): Path to the trajectory XYZ file.
--solv_charge (-c): Path to the solvent charge file.
--qm_radius (-r_qm): QM radius in Å. Solvent within this radius from any dye atom are included in the QM region (default: 5.0 Å).
--mm_radius (-r_mm): MM radius in angstroms (not directly used in the script but can be adjusted as needed).
--nDyes (-n_dyes): Number of dyes in the system.
--dye_atoms (-d_atoms): List of atom counts for each dye. For example, --dye_atoms 17 42 for two dyes with 17 and 42 atoms respectively.
--total_nDyes_atoms (-tot_d_atoms): Total number of atoms in all dyes (sum of atom counts in --dye_atoms).
--nAtoms_solvent (-n_solvent_atoms): Number of atoms per solvent molecule.
--total_frames (-f): Total number of frames (snapshots) in the trajectory.
--total_atoms (-a): Total number of atoms per frame in the trajectory file.
```

### Note

The trajectory file(XYZ) must follow the specified format:

* First line: Total number of atoms per frame
* Second line: Comment or empty.
* Subsequent lines: Atom label and XYZ coordinates, with dyes first and then solvents.
* The solvent charge file (--solv_charge) should contain the charge values for one solvent molecule, with charges listed per atom.
* The script creates a separate directory for each frame, containing the QM and MM files and a TeraChem input file.
* The ordering of atoms in the output files matches the ordering in the input trajectory file.

### Example

Suppose you have:

* A trajectory file named trajectory.xyz with 1000 atoms per frame, including dyes and solvents.
* A solvent charge file named charges.txt containing charges for each atom in a solvent molecule.
* Two dyes in the system: the first dye has 17 atoms, and the second dye has 42 atoms.
* Each solvent molecule consists of 10 atoms.
* You want to process 2 frames

### Output

* The script will generate directories named 1, 2, ..., corresponding to each frame.
* Inside each directory, you will find:
* solute_solvent_5.0ang_qm_gs.xyz: QM region coordinates file.
* solute_solvent_5.0ang_mm.xyz: MM point charges file.
* tc_camb3lyp_5.0ang_opt_gs.in: TeraChem input file.

### Author

Ajay Khanna | Isborn Lab | UC Merced

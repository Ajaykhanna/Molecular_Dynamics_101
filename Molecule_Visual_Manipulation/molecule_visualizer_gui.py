import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import argparse
import os
import sys

# Average bond lengths in Angstroms
average_bond_lengths = {
    frozenset(['C', 'C']): 1.53,
    frozenset(['C', 'N']): 1.47,
    frozenset(['C', 'O']): 1.42,
    frozenset(['C', 'H']): 1.09,
    frozenset(['N', 'H']): 1.00,
    frozenset(['O', 'H']): 0.96,
    frozenset(['C', 'C double']): 1.34,
    frozenset(['C', 'N double']): 1.27,
    frozenset(['C', 'O double']): 1.21,
    frozenset(['C', 'C triple']): 1.20,
    frozenset(['C', 'N triple']): 1.15,
    # Add other bond types as necessary
}

# Tolerance in Angstroms
bond_tolerance = 0.3  # Adjust this value as needed

def get_bond_distance_range(atom1, atom2):
    """
    Retrieve the minimum and maximum bond distances for a given pair of atom types.
    """
    key = frozenset([atom1, atom2])
    average_length = average_bond_lengths.get(key)
    if average_length is not None:
        min_distance = average_length - bond_tolerance
        max_distance = average_length + bond_tolerance
        return min_distance, max_distance
    else:
        return None

def determine_bonds(symbols, coords):
    """
    Determine bonds between atoms based on their types and distances.
    """
    bonds = []
    num_atoms = len(symbols)
    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            atom1 = symbols[i]
            atom2 = symbols[j]
            bond_range = get_bond_distance_range(atom1, atom2)
            if bond_range is not None:
                min_distance, max_distance = bond_range
                distance = np.linalg.norm(coords[i] - coords[j])
                if min_distance <= distance <= max_distance:
                    bonds.append((i, j))
    return bonds

def load_xyz(file_path):
    """
    Load atom symbols and coordinates from an XYZ file.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
        try:
            atom_count = int(lines[0].strip())
        except ValueError:
            raise ValueError("First line must be the number of atoms.")
        data = lines[2:]  # Skip the first two lines (atom count and comment)
        if len(data) < atom_count:
            raise ValueError("Atom count does not match the number of coordinate lines.")
        symbols = []
        coords = []
        for line_num, line in enumerate(data, start=3):
            parts = line.strip().split()
            if not parts:
                continue
            if len(parts) < 4:
                print(f"Warning: Line {line_num} is malformed: '{line.strip()}'")
                continue
            try:
                symbol = parts[0]
                x, y, z = map(float, parts[1:4])
                symbols.append(symbol)
                coords.append([x, y, z])
            except ValueError as e:
                print(f"Error parsing line {line_num}: {e}")
                continue
        return symbols, np.array(coords)

def separate_molecules(symbols, coords, nAtoms_list):
    """
    Separate the atoms and coordinates into molecules based on the number of atoms in each molecule.
    """
    molecules = []
    idx = 0
    for nAtoms in nAtoms_list:
        mol_symbols = symbols[idx:idx + nAtoms]
        mol_coords = coords[idx:idx + nAtoms]
        molecules.append((mol_symbols, mol_coords))
        idx += nAtoms
    return molecules

def calculate_centroid(coords):
    """
    Calculate the centroid of a set of coordinates.
    """
    return np.mean(coords, axis=0)

class MoleculeVisualizer:
    """
    A GUI application for visualizing and manipulating molecules from an XYZ file.
    """
    def __init__(self, master, xyz_file, nMolecules, nAtoms_list):
        """
        Initialize the MoleculeVisualizer.
        """
        self.master = master
        master.title("Molecule Visualizer")
        master.geometry("1000x600")

        self.nMolecules = nMolecules
        self.nAtoms_list = nAtoms_list

        # Allow window resizing
        master.rowconfigure(0, weight=1)
        master.columnconfigure(0, weight=1)

        # Create main frame
        main_frame = tk.Frame(master)
        main_frame.grid(row=0, column=0, sticky="nsew")
        main_frame.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=3)
        main_frame.columnconfigure(1, weight=1)

        # Create plot frame (left side)
        plot_frame = tk.Frame(main_frame)
        plot_frame.grid(row=0, column=0, sticky="nsew")
        plot_frame.rowconfigure(0, weight=1)
        plot_frame.columnconfigure(0, weight=1)

        # Create control frame (right side)
        control_frame = tk.Frame(main_frame)
        control_frame.grid(row=0, column=1, sticky="nsew")
        control_frame.rowconfigure(0, weight=1)
        control_frame.columnconfigure(0, weight=1)

        self.control_frame = control_frame  # Save reference for later use

        # Load molecules
        self.symbols, self.coords = load_xyz(xyz_file)

        # Validate total number of atoms
        if sum(nAtoms_list) != len(self.symbols):
            print("Error: The sum of nAtoms_list does not equal the total number of atoms in the file.")
            sys.exit(1)

        # Separate molecules
        self.molecules = []
        molecule_data = separate_molecules(self.symbols, self.coords, nAtoms_list)
        for mol_symbols, mol_coords in molecule_data:
            mol_bonds = determine_bonds(mol_symbols, mol_coords)
            self.molecules.append({
                'symbols': mol_symbols,
                'coords': mol_coords,
                'bonds': mol_bonds,
                'trans_coords': mol_coords.copy()
            })

        # Default to modifying only molecule 1 (index 0)
        self.selected_molecules = [0]

        # Set up matplotlib figure
        self.fig = plt.figure(figsize=(6, 6))
        self.ax = self.fig.add_subplot(111, projection='3d')

        # Embed the matplotlib figure in Tkinter
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Add controls
        self.add_controls()

        # Initial plot
        self.update_plot()

    def add_controls(self):
        """
        Add translation, rotation, and molecule selection controls to the GUI.
        """
        control_frame = self.control_frame

        # Use a scrollbar if needed
        canvas = tk.Canvas(control_frame)
        scrollbar = tk.Scrollbar(control_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = tk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Molecule Selection
        tk.Label(scrollable_frame, text="Select Molecules to Modify", font=("Arial", 12, "bold")).pack(pady=5)

        # Variable to track "All" selection
        self.select_all_var = tk.IntVar(value=0)  # All unselected by default
        chk_all = tk.Checkbutton(scrollable_frame, text="All Molecules", variable=self.select_all_var, command=self.update_selection)
        chk_all.pack(anchor='w')

        # Variables for individual molecule selections
        self.selection_vars = []
        for i in range(self.nMolecules):
            if i == 0:
                var = tk.IntVar(value=1)  # Select molecule 1 by default
            else:
                var = tk.IntVar(value=0)  # Unselect other molecules
            chk = tk.Checkbutton(scrollable_frame, text=f"Molecule {i+1}", variable=var, command=self.update_selection)
            chk.pack(anchor='w')
            self.selection_vars.append(var)

        # Translation sliders
        self.trans_x = tk.DoubleVar(value=0)
        self.trans_y = tk.DoubleVar(value=0)
        self.trans_z = tk.DoubleVar(value=0)

        tk.Label(scrollable_frame, text="Translation Controls", font=("Arial", 12, "bold")).pack(pady=5)

        tk.Label(scrollable_frame, text="Translate X").pack()
        tk.Scale(scrollable_frame, variable=self.trans_x, orient=tk.HORIZONTAL, length=200, from_=-20, to=20, resolution=0.1, command=self.update_plot).pack()

        tk.Label(scrollable_frame, text="Translate Y").pack()
        tk.Scale(scrollable_frame, variable=self.trans_y, orient=tk.HORIZONTAL, length=200, from_=-20, to=20, resolution=0.1, command=self.update_plot).pack()

        tk.Label(scrollable_frame, text="Translate Z").pack()
        tk.Scale(scrollable_frame, variable=self.trans_z, orient=tk.HORIZONTAL, length=200, from_=-20, to=20, resolution=0.1, command=self.update_plot).pack()

        # Rotation sliders
        self.rot_x = tk.DoubleVar(value=0)
        self.rot_y = tk.DoubleVar(value=0)
        self.rot_z = tk.DoubleVar(value=0)

        tk.Label(scrollable_frame, text="Rotation Controls", font=("Arial", 12, "bold")).pack(pady=5)

        tk.Label(scrollable_frame, text="Rotate X (°)").pack()
        tk.Scale(scrollable_frame, variable=self.rot_x, orient=tk.HORIZONTAL, length=200, from_=-180, to=180, resolution=1, command=self.update_plot).pack()

        tk.Label(scrollable_frame, text="Rotate Y (°)").pack()
        tk.Scale(scrollable_frame, variable=self.rot_y, orient=tk.HORIZONTAL, length=200, from_=-180, to=180, resolution=1, command=self.update_plot).pack()

        tk.Label(scrollable_frame, text="Rotate Z (°)").pack()
        tk.Scale(scrollable_frame, variable=self.rot_z, orient=tk.HORIZONTAL, length=200, from_=-180, to=180, resolution=1, command=self.update_plot).pack()

        # Buttons
        tk.Label(scrollable_frame, text="Actions", font=("Arial", 12, "bold")).pack(pady=5)

        tk.Button(scrollable_frame, text="Save Adjusted Molecule", command=self.save_molecule).pack(pady=5)
        tk.Button(scrollable_frame, text="Reset", command=self.reset_view).pack(pady=5)

    def update_selection(self):
        """
        Update the list of selected molecules based on the user's selection.
        """
        # Update selected molecules list
        self.selected_molecules = [i for i, var in enumerate(self.selection_vars) if var.get() == 1]

        # Update "All Molecules" checkbox based on selections
        if len(self.selected_molecules) == self.nMolecules:
            self.select_all_var.set(1)
        else:
            self.select_all_var.set(0)

        # If "All Molecules" checkbox is toggled
        if self.select_all_var.get():
            for var in self.selection_vars:
                var.set(1)
            self.selected_molecules = list(range(self.nMolecules))

        self.update_plot()

    def update_plot(self, event=None):
        """
        Update the molecule visualization based on current slider values.
        """
        self.ax.clear()

        # Get current transformation parameters
        trans_vector = np.array([self.trans_x.get(), self.trans_y.get(), self.trans_z.get()])
        angle_x = self.rot_x.get()
        angle_y = self.rot_y.get()
        angle_z = self.rot_z.get()

        for idx, molecule in enumerate(self.molecules):
            # Determine if molecule is selected for modification
            if idx in self.selected_molecules:
                # Apply transformation
                coords = molecule['coords']
                # Apply translation
                coords_translated = coords + trans_vector
                # Apply rotation
                coords_transformed = self.apply_rotation(coords_translated, angle_x, angle_y, angle_z)
                molecule['trans_coords'] = coords_transformed
            else:
                # Use original coordinates
                molecule['trans_coords'] = molecule['coords']

            # Plot the molecule
            x, y, z = molecule['trans_coords'].T
            color = plt.cm.tab10(idx % 10)
            self.ax.scatter(x, y, z, color=color, label=f"Molecule {idx+1}")

            # Plot bonds
            for bond in molecule['bonds']:
                i, j = bond
                x_vals = [molecule['trans_coords'][i, 0], molecule['trans_coords'][j, 0]]
                y_vals = [molecule['trans_coords'][i, 1], molecule['trans_coords'][j, 1]]
                z_vals = [molecule['trans_coords'][i, 2], molecule['trans_coords'][j, 2]]
                self.ax.plot(x_vals, y_vals, z_vals, color=color)

        self.ax.legend()
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')
        self.ax.set_title('Molecule Visualization')
        self.ax.view_init(elev=20, azim=30)
        self.canvas.draw()

    def apply_rotation(self, coords, angle_x, angle_y, angle_z):
        """
        Apply rotation to the coordinates based on the provided rotation angles.
        """
        # Convert angles to radians
        theta_x = np.radians(angle_x)
        theta_y = np.radians(angle_y)
        theta_z = np.radians(angle_z)

        # Rotation matrices
        R_x = np.array([[1, 0, 0],
                        [0, np.cos(theta_x), -np.sin(theta_x)],
                        [0, np.sin(theta_x), np.cos(theta_x)]])

        R_y = np.array([[np.cos(theta_y), 0, np.sin(theta_y)],
                        [0, 1, 0],
                        [-np.sin(theta_y), 0, np.cos(theta_y)]])

        R_z = np.array([[np.cos(theta_z), -np.sin(theta_z), 0],
                        [np.sin(theta_z), np.cos(theta_z), 0],
                        [0, 0, 1]])

        # Combined rotation matrix
        R = R_z @ R_y @ R_x

        # Rotate coordinates around the centroid of the molecule
        centroid = calculate_centroid(coords)
        coords_centered = coords - centroid
        coords_rotated = coords_centered @ R.T
        coords_final = coords_rotated + centroid

        return coords_final

    def save_molecule(self):
        """
        Save the adjusted molecule coordinates to an XYZ file.
        """
        symbols = []
        coords = []
        for idx, molecule in enumerate(self.molecules):
            # Use the transformed coordinates
            symbols.extend(molecule['symbols'])
            coords.append(molecule['trans_coords'])
        coords = np.vstack(coords)

        # Write to a new .xyz file
        with open('adjusted_molecule.xyz', 'w') as file:
            file.write(f"{len(symbols)}\n\n")
            for symbol, coord in zip(symbols, coords):
                file.write(f"{symbol}\t{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\n")
        print("Adjusted molecule saved to 'adjusted_molecule.xyz'")

    def reset_view(self):
        """
        Reset the translation and rotation sliders to their initial values.
        """
        # Reset translation sliders to zero
        self.trans_x.set(0)
        self.trans_y.set(0)
        self.trans_z.set(0)

        # Reset rotation sliders to zero
        self.rot_x.set(0)
        self.rot_y.set(0)
        self.rot_z.set(0)

        # Update the plot
        self.update_plot()

def main():
    """
    The main function to parse arguments and start the GUI application.
    """
    parser = argparse.ArgumentParser(description='Visualize and manipulate molecules from an XYZ file.')
    parser.add_argument('xyz_file', type=str, help='Path to the input XYZ file.')
    parser.add_argument('--nMolecules', type=int, required=True, help='Number of molecules in the XYZ file.')
    parser.add_argument('--nAtoms', type=int, nargs='+', required=True, help='Number of atoms in each molecule.')

    args = parser.parse_args()

    # Validate the input file
    if not os.path.isfile(args.xyz_file):
        print(f"Error: File '{args.xyz_file}' does not exist.")
        sys.exit(1)

    if not args.xyz_file.lower().endswith('.xyz'):
        print("Error: Input file must have a .xyz extension.")
        sys.exit(1)

    # Validate nAtoms
    if len(args.nAtoms) != args.nMolecules:
        print("Error: The number of atoms specified does not match the number of molecules.")
        sys.exit(1)

    # Start the GUI application
    root = tk.Tk()
    app = MoleculeVisualizer(root, args.xyz_file, args.nMolecules, args.nAtoms)
    root.mainloop()

if __name__ == "__main__":
    main()

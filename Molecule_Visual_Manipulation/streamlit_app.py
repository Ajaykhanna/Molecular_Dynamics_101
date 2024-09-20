import streamlit as st
import matplotlib.pyplot as plt
import numpy as np

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

def load_xyz_from_text(xyz_text):
    """
    Load atom symbols and coordinates from XYZ file content provided as text.
    """
    lines = xyz_text.strip().splitlines()
    try:
        atom_count = int(lines[0].strip())
    except (ValueError, IndexError):
        st.error("First line must be the number of atoms.")
        return None, None
    data = lines[2:]  # Skip the first two lines (atom count and comment)
    if len(data) < atom_count:
        st.error("Atom count does not match the number of coordinate lines.")
        return None, None
    symbols = []
    coords = []
    for line_num, line in enumerate(data, start=3):
        parts = line.strip().split()
        if not parts:
            continue
        if len(parts) < 4:
            st.warning(f"Line {line_num} is malformed: '{line.strip()}'")
            continue
        try:
            symbol = parts[0]
            x, y, z = map(float, parts[1:4])
            symbols.append(symbol)
            coords.append([x, y, z])
        except ValueError as e:
            st.warning(f"Error parsing line {line_num}: {e}")
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

def apply_rotation(coords, angle_x, angle_y, angle_z):
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

def main():
    st.title("🧪 Molecule Visualizer")

    # Input selection: Upload file or paste content
    input_option = st.radio("Select input method:", ("Upload XYZ File", "Paste XYZ Content"))

    if input_option == "Upload XYZ File":
        uploaded_file = st.file_uploader("Upload an XYZ file", type=["xyz"])
        if uploaded_file is not None:
            try:
                file_content = uploaded_file.getvalue().decode("utf-8")
                symbols, coords = load_xyz_from_text(file_content)
                if symbols is None or coords is None:
                    return
            except Exception as e:
                st.error(f"Error reading file: {e}")
                return
        else:
            st.info("Please upload an XYZ file to proceed.")
            return
    else:
        st.write("Paste your XYZ file content in the text area below:")
        xyz_text = st.text_area("XYZ File Content", height=200)
        if xyz_text:
            symbols, coords = load_xyz_from_text(xyz_text)
            if symbols is None or coords is None:
                return
        else:
            st.info("Please paste the XYZ file content to proceed.")
            return

    # Get number of molecules and atoms
    st.sidebar.header("Molecule Configuration")

    nMolecules = st.sidebar.number_input("Number of Molecules", min_value=1, step=1, value=1)
    nAtoms_input = st.sidebar.text_input("Number of Atoms in Each Molecule (comma-separated)", value="")

    if nAtoms_input:
        try:
            nAtoms_list = [int(x.strip()) for x in nAtoms_input.split(',')]
            if len(nAtoms_list) != nMolecules:
                st.error("Number of atoms specified does not match the number of molecules.")
                return
        except ValueError:
            st.error("Please enter valid integers for the number of atoms.")
            return

        # Validate total number of atoms
        if sum(nAtoms_list) != len(symbols):
            st.error("The sum of atoms does not equal the total number of atoms in the file.")
            return

        # Separate molecules
        molecules = []
        molecule_data = separate_molecules(symbols, coords, nAtoms_list)
        for mol_symbols, mol_coords in molecule_data:
            mol_bonds = determine_bonds(mol_symbols, mol_coords)
            molecules.append({
                'symbols': mol_symbols,
                'coords': mol_coords,
                'bonds': mol_bonds,
                'trans_coords': mol_coords.copy()
            })

        # Molecule Selection
        st.sidebar.header("Select Molecules to Modify")
        selection_options = [f"Molecule {i+1}" for i in range(nMolecules)]
        default_selection = [selection_options[0]]  # Select Molecule 1 by default
        selected_molecules = st.sidebar.multiselect("Molecules", selection_options, default=default_selection)
        selected_indices = [int(s.split()[-1])-1 for s in selected_molecules]

        # Transformation Controls
        st.sidebar.header("Transformation Controls")
        trans_x = st.sidebar.slider("Translate X", -20.0, 20.0, 0.0, 0.1)
        trans_y = st.sidebar.slider("Translate Y", -20.0, 20.0, 0.0, 0.1)
        trans_z = st.sidebar.slider("Translate Z", -20.0, 20.0, 0.0, 0.1)
        rot_x = st.sidebar.slider("Rotate X (°)", -180, 180, 0, 1)
        rot_y = st.sidebar.slider("Rotate Y (°)", -180, 180, 0, 1)
        rot_z = st.sidebar.slider("Rotate Z (°)", -180, 180, 0, 1)

        # Apply transformations and plot
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')

        for idx, molecule in enumerate(molecules):
            # Determine if molecule is selected for modification
            if idx in selected_indices:
                # Apply transformation
                coords = molecule['coords']
                # Apply translation
                coords_translated = coords + np.array([trans_x, trans_y, trans_z])
                # Apply rotation
                coords_transformed = apply_rotation(coords_translated, rot_x, rot_y, rot_z)
                molecule['trans_coords'] = coords_transformed
            else:
                # Use original coordinates
                molecule['trans_coords'] = molecule['coords']

            # Plot the molecule
            x, y, z = molecule['trans_coords'].T
            color = plt.cm.tab10(idx % 10)
            ax.scatter(x, y, z, color=color, label=f"Molecule {idx+1}")

            # Plot bonds
            for bond in molecule['bonds']:
                i, j = bond
                x_vals = [molecule['trans_coords'][i, 0], molecule['trans_coords'][j, 0]]
                y_vals = [molecule['trans_coords'][i, 1], molecule['trans_coords'][j, 1]]
                z_vals = [molecule['trans_coords'][i, 2], molecule['trans_coords'][j, 2]]
                ax.plot(x_vals, y_vals, z_vals, color=color)

        ax.legend()
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Molecule Visualization')
        st.pyplot(fig)

        # Save adjusted molecule
        if st.button("Save Adjusted Molecule"):
            symbols_combined = []
            coords_combined = []
            for molecule in molecules:
                symbols_combined.extend(molecule['symbols'])
                coords_combined.append(molecule['trans_coords'])
            coords_combined = np.vstack(coords_combined)

            # Create XYZ content
            xyz_content = f"{len(symbols_combined)}\n\n"
            for symbol, coord in zip(symbols_combined, coords_combined):
                xyz_content += f"{symbol}\t{coord[0]:.6f}\t{coord[1]:.6f}\t{coord[2]:.6f}\n"

            # Provide download link
            st.download_button(
                label="Download Adjusted Molecule",
                data=xyz_content,
                file_name='adjusted_molecule.xyz',
                mime='text/plain'
            )
    else:
        st.info("Please enter the number of atoms for each molecule in the sidebar.")

if __name__ == "__main__":
    main()

# 🧪 Molecule Visualizer

Welcome to the **Molecule Visualizer**! This interactive GUI application allows you to visualize and manipulate molecules from an XYZ file. With this tool, you can translate and rotate selected molecules, adjust bond criteria, and save the modified configurations.

![Molecule Visualizer Banner](https://your-image-link.com/banner.png)

---

## 🌟 Features

- **Multiple Molecule Support**: Visualize any number of molecules from an XYZ file.
- **Selective Transformation**: Choose which molecules to modify with translation and rotation.
- **Interactive GUI**: User-friendly interface built with Tkinter and Matplotlib.
- **Bond Visualization**: Bonds are drawn based on customizable bond length criteria.
- **Save Configurations**: Export the adjusted molecule structures to a new XYZ file.
- **Resizable Window**: The interface adapts to different screen sizes.

---

## 🖥️ Demo

[![Demo Video](https://img.youtube.com/vi/your-video-id/0.jpg)](https://www.youtube.com/watch?v=your-video-id)

---

## 📦 Installation

### Prerequisites

- **Python 3.x**
- **pip** package manager

### Required Python Packages

Install the required packages using pip:

```bash
pip install numpy matplotlib
```

> **Note**: Tkinter is included with most Python installations. If not, install it via your system's package manager.

### Clone the Repository

```bash
git clone https://github.com/yourusername/molecule-visualizer.git
cd molecule-visualizer
```

---

## 🚀 Usage

### Command-Line Arguments

- `xyz_file`: Path to the input XYZ file.
- `--nMolecules`: Number of molecules in the XYZ file.
- `--nAtoms`: Number of atoms in each molecule (space-separated list).

### Example

```bash
python molecule_visualizer.py molecule.xyz --nMolecules 3 --nAtoms 42 17 29
```

### Running the Application

1. Ensure your XYZ file is correctly formatted.
2. Run the script with the appropriate arguments.
3. The GUI window will appear, displaying the molecules.

---

## 🎛️ GUI Controls

### Molecule Selection

- **Default**: Only **Molecule 1** is selected for modification.
- **Select Molecules to Modify**: Use the checkboxes to select or unselect molecules.
- **All Molecules**: Select or unselect all molecules at once.

### Transformation Controls

- **Translation Controls**:
  - **Translate X**: Move selected molecules along the X-axis.
  - **Translate Y**: Move selected molecules along the Y-axis.
  - **Translate Z**: Move selected molecules along the Z-axis.
- **Rotation Controls**:
  - **Rotate X (°)**: Rotate selected molecules around the X-axis.
  - **Rotate Y (°)**: Rotate selected molecules around the Y-axis.
  - **Rotate Z (°)**: Rotate selected molecules around the Z-axis.

### Actions

- **Save Adjusted Molecule**: Save the current configuration to `adjusted_molecule.xyz`.
- **Reset**: Reset all sliders and selections to their default values.

---

## 📊 Visualization

- **Atoms**: Displayed as colored spheres. Each molecule has a unique color.
- **Bonds**: Lines drawn between atoms based on bond length criteria.
- **Interactive Plot**: Zoom, pan, and rotate the 3D visualization.

---

## ⚙️ Configuration

### Adjusting Bond Tolerance

The `bond_tolerance` variable in the script defines the ± tolerance for bond lengths:

```python
bond_tolerance = 0.3  # Adjust this value as needed
```

- **Decrease**: For stricter bond detection.
- **Increase**: To allow more variation in bond lengths.

### Adding Bond Types

Add more bond types to the `average_bond_lengths` dictionary:

```python
average_bond_lengths = {
    # Existing bond types...
    frozenset(['N', 'O']): 1.40,  # Example average bond length
    # Add other bond types as necessary
}
```

---

## 📝 XYZ File Format

Ensure your XYZ file follows the standard format:

```
<Number of atoms>
<Comment line (can be empty)>
<Atom 1 symbol> <X> <Y> <Z>
<Atom 2 symbol> <X> <Y> <Z>
...
```

- **Atoms**: List all atoms with their symbols and coordinates.
- **Molecules**: Concatenate molecule data sequentially in the file.

---

## ❓ FAQ

### **Q**: I get an error about atom counts not matching.

**A**: Ensure that the sum of atoms specified in `--nAtoms` equals the total number of atoms in your XYZ file.

### **Q**: The GUI doesn't display properly.

**A**: Make sure all dependencies are installed, and you're running Python 3.x.

### **Q**: Bonds aren't displayed correctly.

**A**: Adjust the `bond_tolerance` or verify that the atom types are included in `average_bond_lengths`.

---

## 🤝 Contributing

Contributions are welcome! Please open an issue or submit a pull request.

---

## 📄 License

This project is licensed under the **MIT License**. See the [LICENSE](LICENSE) file for details.

---

## 📬 Contact

- **Author**: Your Name
- **Email**: your.email@example.com
- **GitHub**: [yourusername](https://github.com/yourusername)

---

## ⭐ Acknowledgments

- Thanks to the open-source community for the tools and libraries.
- Inspired by molecular visualization tools in computational chemistry.

---

## 🛠️ Built With

- [Python](https://www.python.org/)
- [Tkinter](https://docs.python.org/3/library/tkinter.html)
- [Matplotlib](https://matplotlib.org/)
- [NumPy](https://numpy.org/)

---

## 🎉 Enjoy Visualizing Your Molecules!
Feel free to reach out if you have any questions or need assistance.
Happy visualizing!


*Made with ❤️ by Ajay Khanna*
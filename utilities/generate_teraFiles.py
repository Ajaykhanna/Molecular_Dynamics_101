import os


def default_tc_params(precision="double", dft_grid=5):
    """
    TeraChem Default input parameters.
    """

    return f"""# Job Accuracy & Hardware Control
precision     {precision}
threall       1.0e-15
convthre      3.0e-7
dftgrid       {dft_grid}
gpus          all
safemode      no
"""


def default_tc_excited_state_params(root, nstates, nvecs, density_difference):
    """
    TeraChem Default input parameters for excited state calculations.
    """

    return f"""# Excited State Parameters
cis            yes
cistarget      {root}
cisnumstates   {nstates}
cisguessvecs   {nvecs}
cisdiffdensity {density_difference}
cismaxiter     55
"""


def get_constaint_freeze(nDyes_atoms, coordinate_file, opt_constraints=True):
    """
    Generate TeraChem constraint/freeze input parameters based on the QM region.
    
    Args:
        nDyes_atoms (int): The number of atoms in the dye region.
        coordinate_file (str): The filename of the QM coordinates file.
        opt_constraints (bool): Whether to generate constraints or freeze the non-dye region.
    
    Returns:
        str: A string containing the TeraChem constraint/freeze input parameters.
    """
    newline = "\n"
    with open(coordinate_file, "r") as f:
        lines = f.readlines()
        total_qm_atoms = int(lines[0])
    
    if opt_constraints:
        return f"""# Freezing Atom's XYZ Coordinates
$constraint_freeze
xyz {nDyes_atoms + 1}-{total_qm_atoms}
$end
"""
    else:
        return f"""# Freezing Atoms
$constraints
{''.join(f'atom {i}{newline}' for i in range(nDyes_atoms + 1, total_qm_atoms + 1))}
$end
"""


def generate_tc_vertical_excitation_energy_file(
    point_charge_file,
    coordinate_file,
    net_charge=1,
    spin_mult=1,
    XC_functional="camb3lyp",
    basis_set="6-31G*",
    root=1,
    nstates=6,
    nvecs=6,
    density_difference="yes",
    scr_dir="vee",
):
    """
    Generate TeraChem input parameters for a vertical excitation energy calculation.
    
    Args:
        point_charge_file (str): Filename of the MM point charges file.
        coordinate_file (str): Filename of the QM coordinates file.
        net_charge (int): The net charge of the system.
        spin_mult (int): The spin multiplicity of the system.
        XC_functional (str): The exchange-correlation functional to use.
        basis_set (str): The basis set to use.
        root (int): The root state to target for the excitation.
        nstates (int): The number of excited states to calculate.
        nvecs (int): The number of guess vectors to use.
        density_difference (str): Whether to calculate the density difference.
        scr_dir (str): The scratch directory to use for the calculation.
    
    Returns:
        str: A string containing the TeraChem input parameters.
    """
    return f"""# TeraChem Job-Control Info.
charge         {net_charge}
spinmult       {spin_mult}
basis          {basis_set}
method         {XC_functional}

# Type of Job: Vertical Excitation Energy
run           energy

{default_tc_excited_state_params(root, nstates, nvecs, density_difference)}

# XYZ Filename/path
pointcharges  {os.path.basename(point_charge_file)}
coordinates   {os.path.basename(coordinate_file)}
dispersion    d3

# Scratch Directory Information
scrdir        {scr_dir}
keep_scr      yes

{default_tc_params()}
end
"""


def generate_tc_ground_state_optimization_file(
    point_charge_file,
    coordinate_file,
    nDyes_atoms=0,
    net_charge=0,
    spin_mult=1,
    XC_functional="camb3lyp",
    basis_set="6-31G*",
):
    """
Generate a TeraChem input file for a ground state optimization calculation.

Args:
    point_charge_file (str): Filename of the MM point charges file.
    coordinate_file (str): Filename of the QM coordinates file.
    nDyes_atoms (int): The number of atoms in the dye molecules.
    net_charge (int): The net charge of the system.
    spin_mult (int): The spin multiplicity of the system.
    XC_functional (str): The exchange-correlation functional to use.
    basis_set (str): The basis set to use.

Returns:
    str: A string containing the TeraChem input parameters for a ground state optimization.
"""

    return f"""# TeraChem Job-Control Info.
charge         {net_charge}
spinmult       {spin_mult}
basis          {basis_set}
method         {XC_functional}

# Type of Job: Ground State Optimization
run           minimize
new_minimizer yes

# XYZ Filename/path
pointcharges  {os.path.basename(point_charge_file)}
coordinates   {os.path.basename(coordinate_file)}
dispersion    d3

# Scratch Directory Information
scrdir        gs_opt
keep_scr      yes

{default_tc_params()}
end

{get_constaint_freeze(nDyes_atoms, coordinate_file, opt_constraints=True)}
"""

def generate_tc_ground_state_frequency_file(
    point_charge_file,
    coordinate_file,
    nDyes_atoms=0,
    net_charge=0,
    spin_mult=1,
    XC_functional="camb3lyp",
    basis_set="6-31G*",
):
    """
Generate a TeraChem input file for a ground state optimization calculation.

Args:
    point_charge_file (str): Filename of the MM point charges file.
    coordinate_file (str): Filename of the QM coordinates file.
    nDyes_atoms (int): The number of atoms in the dye molecules.
    net_charge (int): The net charge of the system.
    spin_mult (int): The spin multiplicity of the system.
    XC_functional (str): The exchange-correlation functional to use.
    basis_set (str): The basis set to use.

Returns:
    str: A string containing the TeraChem input parameters for a ground state optimization.
"""

    return f"""# TeraChem Job-Control Info.
charge         {net_charge}
spinmult       {spin_mult}
basis          {basis_set}
method         {XC_functional}

# Type of Job: Ground State Optimization
run           frequencies
mincheck      false

# XYZ Filename/path
pointcharges  {os.path.basename(point_charge_file)}
coordinates   {os.path.basename(coordinate_file)}
dispersion    d3

# Scratch Directory Information
scrdir        gs_freq
keep_scr      yes

{default_tc_params()}
end

{get_constaint_freeze(nDyes_atoms, coordinate_file, opt_constraints=False)}
"""


if __name__ == "__main__":
    print("Script to Generate TeraChem Input Files")

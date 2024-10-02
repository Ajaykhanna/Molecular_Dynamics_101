import os

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
charge         0
spinmult       1
basis          6-31G*
method         camb3lyp

# Excited State Calculation
cis            yes
cistarget      1
cisnumstates   6
cisguessvecs   6
cisdiffdensity yes
cismaxiter     55

# Type of Job: Vertical Excitation Energy
run           energy

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
    print("Script to Generate TeraChem Input Files")
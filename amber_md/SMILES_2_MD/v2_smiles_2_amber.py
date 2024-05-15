"""
High-throughput Molecular Dynamics Simulations From SMILES
of drug-like small molecules for virtual screening.

Author: Ajay Khanna
Date: Jun.15.2023
Place: Frontier Medicines, SSF
"""

import os
import argparse
import pandas as pd
from tqdm import tqdm
import subprocess as subp
from multiprocessing import Pool, cpu_count

def smiles_2_pdb(args):
    """
    Converts a SMILES string to a PDB file using OpenBabel and pdb4amber.

    Args:
        args (tuple): A tuple containing the SMILES string and the output file name.

    Returns:
        None
    """
    smiles, ofile = args
    unprocessed_pdb = subp.getoutput(
        f"obabel -:'{smiles}' -O {ofile}/unprocessed_drug.pdb --gen3d"
    )
    processed_pdb = subp.getoutput(
        f"pdb4amber -i {ofile}/unprocessed_drug.pdb -o {ofile}/{ofile}.pdb"
    )

    return None


def pdb_2_prmtop(args):
    """
    Converts a PDB file to Amber parameter files (prmtop and inpcrd) using the tleap program.

    Args:
        args (str): The name of the PDB file to be converted.

    Returns:
        None
    """
    ofile = args
    os.chdir(f"{ofile}")
    pdb_2_mol2 = subp.getoutput(
        f"antechamber -i {ofile}.pdb -fi pdb -o {ofile}.mol2 -fo mol2 -c bcc -s 2 -at gaff -nc 0 -pl 30"
    )
    mol2_parmchk = subp.getoutput(f"parmchk2 -i {ofile}.mol2 -f mol2 -o {ofile}.frcmod")
    tleap_prmtop = subp.getoutput(
        f"""
    cat > tleap_input.inp << EOF
# Loading FFs: TIP3P and GAFF
source leaprc.water.tip3p
source leaprc.gaff

# Loading Ligands Amber's Parameters
mol = loadmol2 {ofile}.mol2
loadamberparams {ofile}.frcmod
check mol

# Saving Ligand's Parameters: Gas Phase
mol_gas = copy mol
saveamberparm mol_gas ligand_gas.prmtop ligand_gas.inpcrd
savepdb mol_gas ligand_gas.pdb

# Saving Ligand's Parameters: Solvent Phase
mol_water = copy mol

# Setting an Initial Size of the Box
# This Will Change Latter on to Adjust for Water
# Molecules in a PBC ~Cubic Box
set mol_water box 15.0

# Setting Box Center as Ligand's Geometric Center
setbox mol_water centers

# Solvating the Ligand Box with pre-equilibrated TIP3PBOX
# 15.870 = Distance of the Ligand Center from the Edges in all direction
# This is also called the buffer size
# Iso: x = y = z = cube
solvatebox mol_water TIP3PBOX 15.870 iso
saveamberparm mol_water ligand_water.prmtop ligand_water.inpcrd
savepdb mol_water ligand_water.pdb

quit
EOF
"""
    )
    tleap_run = subp.getoutput(f"tleap -f tleap_input.inp")
    os.chdir("../")

    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Python program to convert SMILES to Amber Parameter Files"
    )
    parser.add_argument(
        "--ifile", type=str, required=True, help="SMILES file in CSV format"
    )
    parser.add_argument("--ofile", type=str, required=True, help="Output pdb root name")
    parser.add_argument(
        "--row1",
        type=int,
        default=0,
        required=False,
        help="Number of molecules to process, default=First molecule",
    )
    parser.add_argument(
        "--row2",
        type=int,
        default=1,
        required=False,
        help="Number of molecules to process, default=1 molecule",
    )
    args = parser.parse_args()

    drugs_df = pd.read_csv(args.ifile, header=0)
    drug_name = [
        f"{args.ofile}_drug_{index}"
        for index in range(len(drugs_df))[args.row1 : args.row2]
    ]
    smiles = [smile for smile in drugs_df["smiles"][args.row1 : args.row2]]

    if len(drug_name) == len(smiles):
        with Pool(cpu_count()) as p:
            list(
                tqdm(p.imap(smiles_2_pdb, zip(smiles, drug_name)), total=len(drug_name))
            )
            list(tqdm(p.imap(pdb_2_prmtop, drug_name), total=len(drug_name)))
    else:
        print("Mismatch between #SMILES and # Compound names")

    print("Done")


# How to Run? : python v2_smiles_2_amber.py --ifile smiles.csv --ofile 5VFI --row1 11 --row2 21

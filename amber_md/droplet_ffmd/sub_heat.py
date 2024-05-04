#---------------------------
# Developer: Ajay Khanna
# Place: UC Merced
# Date: April.30.2024
#---------------------------

import sys
import os

infile = sys.argv[1]
name = os.path.splitext(infile)[0]
outfile = name + ".sub"

def generate_minimization_banner():
    """
    Generates a banner for submitting a minimization job using Amber 20.

    Returns:
        str: The banner text.
    """
    banner = """
    +--------------------------------------------------+
    |  Submitting NVT Quick Heating Job Using Amber 20 |
    +--------------------------------------------------+
    """
    return banner


# Open the output file for writing
with open(outfile, 'w') as fileobj:
    fileobj.write("#!/bin/sh\n")
    fileobj.write("#SBATCH --nodes=1\n")
    fileobj.write("#SBATCH --ntasks=32\n")
    fileobj.write("#SBATCH --partition short\n")
    fileobj.write("##SBATCH --mem=200G\n")
    fileobj.write("#SBATCH --time=06:00:00 # 06 hours\n")
    fileobj.write("#SBATCH --output=" + name + ".out\n")
    fileobj.write("#SBATCH --job-name=" + name + "\n")
    fileobj.write("#SBATCH --export=ALL\n")
    fileobj.write("##SBATCH --exclusive\n\n")
    fileobj.write("whoami\n")
    fileobj.write("module load openmpi/4.0.6-gcc-8.4.1\n")
    fileobj.write("module load cuda/11.0.3\n")
    fileobj.write("export AMBERHOME=/home/akhanna2/data/test_amber_modifications/amber20\n")
    fileobj.write("export CUDA_VISIBLE_DEVICES=0,1\n")
    fileobj.write("export CUDA_HOME=/opt/spack/opt/spack/linux-rhel8-icelake/gcc-8.4.1/cuda-11.0.3-mwxjoa2nfce32sfsxxxz66oza5p5tr35/\n")
    fileobj.write("echo $LD_LIBRARY_PATH\n")
    fileobj.write("echo $PATH\n")
    fileobj.write("mpirun -np 32 $AMBERHOME/bin/sander.MPI  -O -i heat.in -o heat.out -p solute_solvent_drop.prmtop -c min.ncrst -r heat.ncrst -x heat.mdcrd -inf heat.info -ref min.ncrst \n")

# Execute the sbatch command
cmd = "sbatch " + outfile
os.system(cmd)
print(generate_minimization_banner())

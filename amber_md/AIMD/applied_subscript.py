#---------------------------
# Developer: Ajay Khanna
# Place: UC Merced
# Date: April.30.2024
#---------------------------

import sys
import os

outfile="mdin.sub"

def generate_aimd_banner():
    """
    Generates a banner for submitting a minimization job using Amber 20.

    Returns:
        str: The banner text.
    """
    banner = """
    +--------------------------------------------------+
    |  Submitting NVT AIMD Job Using Amber 20 	       |
    +--------------------------------------------------+
    """
    return banner

# Open the output file for writing
with open(outfile, 'w') as fileobj:
	fileobj.write("#!/bin/csh\n")
	fileobj.write("#SBATCH --nodes=1\n")
	fileobj.write("#SBATCH --ntasks=32\n")
	fileobj.write("#SBATCH -p test\n")
	fileobj.write("##SBATCH --mem=100G\n")
	fileobj.write("#SBATCH --time=00-01:00:00     # 8 hours\n")
	fileobj.write("#SBATCH --output="+name+".out\n")
	fileobj.write("#SBATCH --job-name="+name+"\n")
	fileobj.write("#SBATCH --export=ALL\n")
	fileobj.write("#SBATCH --gres gpu:1\n")
    fileobj.write("whoami\n")
    fileobj.write("module load gaussian/g16-b01\n")
    fileobj.write("module list\n")
    fileobj.write("echo $LD_LIBRARY_PATH\n")
    fileobj.write("echo $PATH\n")
    fileobj.write("export AMBERHOME=/home/akhanna2/data/amber20\n")
    fileobj.write("./sub_qmmm_terachem.sh\n")
    fileobj.close()
    
# Execute the sbatch command
cmd = "sbatch " + outfile
os.system(cmd)
print(generate_aimd_banner())

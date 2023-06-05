#!/bin/sh
## --------------------------------------------------------
## Information
## @email: akhanna2@ucmerced.edu || quantphobia@gmail.com
## @lab: Dr. Isborn
## @place: UC Merced
## @date: July.10.2021
## @author: Ajay Khanna
## --------------------------------------------------------

# submit_array.sh

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --partition gpu
#SBATCH --mem=10G
#SBATCH --time=03-00:00:00     # 03 Days
#SBATCH --output=openmm_log.out
#SBATCH --open-mode=append     # Append to Outfile
#SBATCH --job-name=OpenMMMD
##SBATCH --exclude=gnode003
#SBATCH --export=ALL
#SBATCH --gres gpu:2

whoami
module load cuda/10.2.89
module list
conda activate openmm_env
export CUDA_VISIBLE_DEVICES=0,1

# Work around slurm arraying
nvidia-smi --query-gpu=index,memory.used --format=csv --loop=4.7 >> OpenMM_Job_GPUsage.data &
python full_openmm.py >> output.log 


echo "All Done"



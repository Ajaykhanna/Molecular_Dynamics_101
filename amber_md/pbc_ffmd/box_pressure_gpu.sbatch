#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH -p test
#SBATCH --gres=gpu:2
#SBATCH --time=00-01:00:00
#SBATCH --output=./job.out
#SBATCH --job-name=amber_pressure
#SBATCH --export=ALL

module load cuda/11.0.3
module load openmpi/4.0.6-gcc-8.4.1
export CUDA_VISIBLE_DEVICES=0,1
export CUDA_HOME=/opt/spack/opt/spack/linux-rhel8-icelake/gcc-8.4.1/cuda-11.0.3-mwxjoa2nfce32sfsxxxz66oza5p5tr35/
export AMBERHOME=/home/akhanna2/data/test_amber_modifications/amber20

whoami
module list
echo $LD_LIBRARY_PATH
echo $PATH

$AMBERHOME/bin/pmemd.cuda -O -i box_pressure.in -o pressure.out -p solute_solvent_box.prmtop  -c heat.ncrst -r pressure.ncrst -x pressure.mdcrd -ref heat.ncrst -inf pressure.info

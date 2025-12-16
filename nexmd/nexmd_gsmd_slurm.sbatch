#!/bin/bash

# Submit this script with: sbatch filename

#SBATCH --time=02-00:00:00   # walltime
#SBATCH --qos=long
#SBATCH --mem=0              # Request all memory on the node if mem=0, use a specific value like --mem=128GB instead
#SBATCH --nodes=1            # number of nodes
#SBATCH --ntasks=4           # number of processor cores (i.e. tasks)
#SBATCH --job-name=ACN_MD    # job name
##SBATCH --account=s17_cint   # account name (uncomment and set your account)
#SBATCH --partition=ml4chem  # partition name
#SBATCH --no-requeue         # do not requeue when preempted and on node failure
#SBATCH --signal=23@60       # send signal to job at [seconds] before end

# Standard production practices
#SBATCH --output=slurm-%j.out   # Standard output and error log
#SBATCH --error=slurm-%j.err    # error log file
#SBATCH --mail-type=END,FAIL    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=user@example.com # Where to send mail (replace with your email)

# LOAD MODULEFILES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

# --- Start Time Tracking ---
START_TIME=$(date +"%Y-%m-%d %H:%M:%S")
echo "Job started on: $START_TIME"
echo "SLURM Job ID: $SLURM_JOB_ID"

# Set robust error handling
set -e
trap 'echo "Job failed with exit code $?"; echo "Job ended on: $(date +"%Y-%m-%d %H:%M:%S")"; exit 1' ERR

ulimit -s unlimited

# Run your program
echo "Running NEXMD simulation..."
/akhanna2/softwares/NEXMD/nexmd.exe < input.ceon > output_gs_md.out

# Check the exit status of the simulation command (if the trap didn't catch it)
if [ $? -eq 0 ]; then
    echo "Simulation completed successfully."
else
    echo "Simulation finished with an error."
fi

sleep 5 # The original script had a sleep command, keeping it here.

# --- End Time Tracking ---
END_TIME=$(date +"%Y-%m-%d %H:%M:%S")
echo "Job ended on: $END_TIME"

# Calculate and display total time consumption
# This uses date command features for time difference calculation
# The exact syntax for date arithmetic can vary by system (e.g., GNU date vs BSD date)
# This example uses a common method for elapsed seconds and simple formatting:
if command -v gdate >/dev/null 2>&1; then
    # Use gdate if available (e.g., on macOS)
    START_SECS=$(gdate -d "$START_TIME" +%s)
    END_SECS=$(gdate -d "$END_TIME" +%s)
else
    # Use standard date command (common on Linux HPCs)
    START_SECS=$(date -d "$START_TIME" +%s)
    END_SECS=$(date -d "$END_TIME" +%s)
fi

ELAPSED_SECS=$((END_SECS - START_SECS))
ELAPSED_TIME=$(printf '%dd-%dh:%dm:%ds\n' $((ELAPSED_SECS/86400)) $((ELAPSED_SECS%86400/3600)) $((ELAPSED_SECS%3600/60)) $((ELAPSED_SECS%60)))

echo "Total time consumed: $ELAPSED_TIME"

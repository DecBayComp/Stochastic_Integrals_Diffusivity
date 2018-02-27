#! /bin/bash

# Sbatch options
#SBATCH -J ito-guy
#SBATCH --qos=fast
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=100MB

# Constants
logs_folder="./logs/"
args_file="arguments.dat"

# Read command line arguments from file
argument=`awk "NR==${SLURM_ARRAY_TASK_ID}" $args_file`

# Launch srun with these argument sequence
echo $argument
srun -o "${logs_folder}log_job_${SLURM_ARRAY_TASK_ID}.out" -e "${logs_folder}log_job_${SLURM_ARRAY_TASK_ID}.err" -J "${SLURM_ARRAY_TASK_ID}" python3 simulate_one_trajectory.py $argument




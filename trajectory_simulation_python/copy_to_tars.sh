#! /bin/bash

FILE_LIST=("batch_tars_simulate_trajectories.sh constants.py D_func.py f_func.py sbatch_one_job.sh simulate_one_trajectory.py")

for file in $FILE_LIST
do
	scp ./$file aserov@tars:~/ito-stratonovich/
done

echo "All files copied successfully"
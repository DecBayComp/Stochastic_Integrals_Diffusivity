#! /bin/bash

FILE_LIST=("batch_t_bayes_simulate_trajectories.sh constants.py D_func.py f_func.py sbatch_one_job_t_bayes.sh simulate_one_trajectory.py")

for file in $FILE_LIST
do
	scp ./$file aserov@thomas_bayes:~/ito-stratonovich/
done

echo "All files copied successfully"
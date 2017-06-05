#! /bin/bash


## Constants
min_D_case_number=3
max_D_case_number=4

min_f_case_number=1
max_f_case_number=5


echo "Submitting tasks to kernels..."

for D_case_number in $(seq $min_D_case_number $max_D_case_number)
do	
	for f_case_number in $(seq $min_f_case_number $max_f_case_number)
	do
		echo "Submitting task with D=${D_case_number}, f=${f_case_number} to kernels..."
		# nohup python2.7 simulate_one_trajectory.py $D_case_number $f_case_number > "log_D=${D_case_number}_f=${f_case_number}.log" &
		# module load Python/2.7.11
		srun -o "log_D=${D_case_number}_f=${f_case_number}.out" -e "log_D=${D_case_number}_f=${f_case_number}.err" -J "D=${D_case_number}_f=${f_case_number}" --cpus-per-task=1 --mem=32GB python2.7 simulate_one_trajectory.py $D_case_number $f_case_number &
		echo "Task submitted"
		sleep 1s
	done
done

echo "All tasks submitted successfully!"




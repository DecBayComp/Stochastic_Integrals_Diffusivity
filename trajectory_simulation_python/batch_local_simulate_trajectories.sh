#! /bin/bash


## Constants
trials=10
sleep_time=0.0
logs_folder="./logs/"
output_folder="./output/"
args_file="./arguments.dat"
D_case=7
f_case=2

echo "Creating arguments list..."

# # Clear the arguments file
# rm $args_file


# Create the logs folder
echo "Cleaning up the logs folder"
if [ ! -d "$logs_folder" ]
then
	mkdir $logs_folder
# Else empty the folder
else
	rm -v ${logs_folder}*
fi


# Create the output folder
echo "Cleaning up the output folder"
if [ ! -d "$output_folder" ]
then
	mkdir $output_folder
# Else empty the folder
else
	rm -v ${output_folder}*
fi



id=0
for ((trial=1;trial<=trials;trial++))
do	
	# Fixed lambda
	for lambda in 0.0 0.5 1.0
	do
		id=$((id+1))
		# echo "Submitting task to cluster. Trial: ${trial}. Lambda: ${lambda}..."
		# module load Python/2.7.11
		# srun -o "${logs_folder}log_lambda_fixed_${id}.out" -e "${logs_folder}log_lambda_fixed_${id}.err" -J "T=${trial}_${lambda}" --cpus-per-task=1 --mem=100MB --qos=fast python simulate_one_trajectory.py -D=6 -f=7 -l=$lambda --id=$id &
		
		`python2 simulate_one_trajectory.py -D=${D_case} -f=${f_case} -l=$lambda --id=$id >> ${logs_folder}log_lambda_fixed_${id}.out &`
		# echo "-D=${D_case} -f=${f_case} -l=$lambda --id=$id" >> $args_file

		
	done

	# Random lambda
	id=$((id+1))
	# echo "Submitting task to cluster. Trial: ${trial}. Lambda: rand..."
	# module load Python/2.7.11
	# srun -o "${logs_folder}log_lambda_rand_${trial}.out" -e "${logs_folder}log_lambda_rand_${trial}.err" -J "T=${trial}_rnd" --cpus-per-task=1 --mem=100MB --qos=fast python simulate_one_trajectory.py -D=6 -f=7 --rand --id=$trial &
	
	# echo "-D=${D_case} -f=${f_case} --rand --id=$id" >> $args_file
	`python2 simulate_one_trajectory.py -D=${D_case} -f=${f_case} --rand --id=$id >> ${logs_folder}log_lambda_rand_${trial}.out &`
	
	# sleep ${sleep_time}s
	# echo "Task submitted"
done
echo "All tasks submitted..."

# lines_count=$id

# # Launch sbatch
# sbatch -o /dev/null --array=1-$lines_count sbatch_one_job.sh


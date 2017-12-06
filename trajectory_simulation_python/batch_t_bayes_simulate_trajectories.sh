#! /bin/bash


## Constants
trials=10
sleep_time=0.2
logs_folder="./logs/"
output_folder="./output/"
args_file="./arguments.dat"
D_case=7
f_force_case=2
f_no_force_case=1

echo "Creating arguments list..."

# Clear the arguments file
rm $args_file

# Create the logs folder
if [ ! -d "$logs_folder" ]
then
	mkdir $logs_folder
# Else empty the folder
else
	rm -v ${logs_folder}*
fi


# Create the output folder
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
		echo "-D=${D_case} -f=${f_force_case} -l=$lambda --id=$id" >> $args_file

		id=$((id+1))	
		echo "-D=${D_case} -f=${f_no_force_case} -l=$lambda --id=$id" >> $args_file
	done

	# Random lambda
	id=$((id+1))
	echo "-D=${D_case} -f=${f_force_case} --rand --id=$id" >> $args_file

	id=$((id+1))
	echo "-D=${D_case} -f=${f_no_force_case} --rand --id=$id" >> $args_file
done
echo "Argument list created. Launching sbatch..."

lines_count=$id

# Launch sbatch
# echo "sbatch -o /dev/null --array=1-$lines_count sbatch_one_job_t_bayes.sh"
sbatch -o /dev/null --array=1-$lines_count sbatch_one_job_t_bayes.sh


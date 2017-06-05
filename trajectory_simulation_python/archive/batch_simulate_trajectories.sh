#! /bin/bash


## Constants
max_D_case_number=4
max_f_case_number=5

for D_case_number in $(seq 1 $max_D_case_number)
do	
	for f_case_number in $(seq 1 $max_f_case_number)
	do
		nohup python simulate_one_trajectory.py $D_case_number $f_case_number > "log_D=${D_case_number}_f=${f_case_number}.log" &
	done
done






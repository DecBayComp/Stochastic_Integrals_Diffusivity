#! /bin/bash





for D_case_number in $(seq 1 5)
do
	srun -J "test_queue" --cpus-per-task=1 --mem=1MB sleep 100 &
done
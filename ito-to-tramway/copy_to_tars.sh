#! /bin/bash

FILE_LIST=("prepare_arguments.py filelock.py LICENSE-filelock.rst start_me.py 
	sbatch_tars.sh job_manager.py tesselate_and_infer.py calculate.py")

for file in $FILE_LIST
do
	scp ./$file aserov@tars.pasteur.fr:~/bayes_factor_illustration/
done

echo "All files copied successfully"
echo

read -p "Press any key to continue... " -n1 -s

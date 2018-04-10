#! /bin/bash

FILE_LIST=("start_me.py constants.py D_func.py main.py sbatch_tars.sh job_manager.py filelock.py")

for file in $FILE_LIST
do
	scp ./$file aserov@tars.pasteur.fr:~/ito-stratonovich/
done

echo "All files copied successfully"
echo

read -p "Press any key to continue... " -n1 -s
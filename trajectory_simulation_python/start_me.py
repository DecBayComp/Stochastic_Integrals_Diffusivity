# This function creates the argument list and launches job managers



import os    # for file operations
import socket	# for netowrk hostname
import numpy as np
import argparse	# for command-line arguments

## Constants
from constants import *


## Define arguments
arg_parser = argparse.ArgumentParser(description = 'Job manager. You must choose whether to resume simulations or restart and regenerate the arguments file')
arg_parser.add_argument('-v', '--version', action = 'version', version = '%(prog)s ' + str(version))
mode_group = arg_parser.add_mutually_exclusive_group(required = True)
mode_group.add_argument('--restart', action = 'store_true')
mode_group.add_argument('--resume', action = 'store_true')


# Identify the system where the code is running
hostname = socket.gethostname()
if hostname.startswith('tars-submit'):
	script_name = 'sbatch_tars.sh'
	jobs_count = jobs_count_tars
elif hostname == 'patmos':
	script_name = 'sbatch_t_bayes.sh'
	jobs_count = jobs_count_t_bayes
else:
	print('Unidentified hostname "' + hostname + '". Unable to choose the right code version to launch. Aborting.')
	exit()


## Analyze if need to restart or resume
input_args = arg_parser.parse_args()
bl_restart = input_args.restart


# If restart is required, regenerate all files
if bl_restart:
	print("Creating arguments list...")

	# Clear the arguments file
	try:
		os.remove(args_file)
	except:
		pass


	# Clean the logs and output folders
	for folder in (logs_folder, output_folder):
		if os.path.isdir(folder):
			print("Cleaning up the folder: '" + folder + "'.")
			cmd = "rm -rfv " + folder 
			try:
				os.system(cmd)
			except Exception as e:
				print(e)

		# Recreate the folder
		try:
			os.makedirs(folder)
		except Exception as e:
			print(e)

	# Clean slurm files in the root folder
	cmd = "rm -fv ./slurm-*" 
	try:
		os.system(cmd)
	except Exception as e:
		print(e)


	## Write arguments to file
	id = 0
	ksi_length = (ksi_range[1] - ksi_range[0]) / ksi_step + 1
	ksi_points = np.linspace(ksi_range[0], ksi_range[1], ksi_length)
	with open(args_file, 'w') as file_object:
		for trial in range(1, trials+1):
			for ksi in ksi_points:
				id += 1
				args_string = '-D=%i --ksi=%f --id=%i\n' % (D_case, ksi, id)
				file_object.write(args_string)


	# Create lock file
	with open(args_lock, 'w'):
		pass


	print("Argument list created. Launching sbatch...")

	line_count = id
else:
	print("Resuming simulation with the same arguments file")


# Launch job managers
cmd_str = 'sbatch --array=1-%i %s' % (jobs_count, script_name)	# -o /dev/null
os.system(cmd_str)


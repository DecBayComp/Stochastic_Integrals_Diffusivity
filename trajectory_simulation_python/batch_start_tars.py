


import os    # for file operations

## Constants
from constants import trials, sleep_time, logs_folder, output_folder, args_file
D_case = 2

print("Creating arguments list...")

# Clear the arguments file
os.remove(args_file)


# Prepare the logs and output folders
for folder in (logs_folder, output_folder):
	print("Cleaning up the folder: " + folder)
	if not os.path.isdir(folder):
		os.makedirs(folder)
		
	# Else empty the folder
	else:
		for file in os.listdir(folder):
			file_path = os.path.join(folder, file)
			try:
				if os.path.isfile(file_path):
					os.remove(file_path)
			except Exception as e:
				print(e)


## Write arguments to file
id = 0
with open(args_file, 'w') as file_object:
	for ksi in range(ksi_range(0), ksi_range(1), ksi_step):
		for trial in range(1, trials+1):
			id += 1
			args_string = '-D=%i --ksi=%f' % (D_case, ksi)
			file_object.write(args_string)

print("Argument list created. Launching sbatch...")

line_count = id

# Launch sbatch
cmd_str = 'sbatch -o /dev/null --array=1-%i sbatch_one_job.sh' % (line_count)
os.system(cmd_str)


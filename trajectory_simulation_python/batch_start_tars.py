


import os    # for file operations
import numpy as np

## Constants
from constants import trials, sleep_time, logs_folder, output_folder, args_file, ksi_range, ksi_step
D_case = 2

print("Creating arguments list...")

# Clear the arguments file
try:
	os.remove(args_file)
except:
	pass


# Clean the logs and output folders
for folder in (logs_folder, output_folder):
	if not os.path.isdir(folder):
		os.makedirs(folder)
		
	# Else empty the folder
	else:
		print("Cleaning up the folder: '" + folder + "'.")
		file_count = len(os.listdir(folder))
		file_num = 0
		for file in os.listdir(folder):
			file_path = os.path.join(folder, file)
			try:
				if os.path.isfile(file_path):
					os.remove(file_path)
					file_num += 1
					print ("File: %i/%i" % (file_num, file_count))
			except Exception as e:
				print(e)


## Write arguments to file
id = 0
ksi_length = (ksi_range[1] - ksi_range[0]) / ksi_step + 1
ksi_points = np.linspace(ksi_range[0], ksi_range[1], ksi_length)
with open(args_file, 'w') as file_object:
	for ksi in ksi_points:
		for trial in range(1, trials+1):
			id += 1
			args_string = '-D=%i --ksi=%f --id=%i\n' % (D_case, ksi, id)
			file_object.write(args_string)

print("Argument list created. Launching sbatch...")

line_count = id

# Launch sbatch
cmd_str = 'sbatch -o /dev/null --array=1-%i sbatch_one_job.sh' % (line_count)
os.system(cmd_str)


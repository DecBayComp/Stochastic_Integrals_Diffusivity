# The function reads one line of arguments from a file and launches code execution.
# It correctly handles the existence of several job_managers in parallel.


from filelock import FileLock, Timeout		# for file locks
import os	# for shell execution of the code
import sys
import socket	# for netowrk hostname
from main import main	# actual calculations

from constants import *
# calculation_script = "main.py"


while True:
	# Try to obtain lock
	lock = FileLock(args_lock, timeout = lock_timeout)
	try:
		# If get lock
		with lock:
			# Read and store all lines in the arguments file
			with open(args_file, 'r') as file_object:
				all_lines = file_object.read().splitlines(True)

			# If arguments file is empty, end program
			if len(all_lines) == 0:
				print("Arguments file empty. Finishing.")
				sys.exit(0)
			
			# Write all lines back except for the first one
			with open(args_file, 'w') as file_object:
				file_object.writelines(all_lines[1:])

			# Get current arguments and clen
			cur_args = all_lines[0].strip()
			del all_lines

	# If unable to get lock
	except Exception as e:
		print("Ecountered unknown exception while getting the file lock: '"+ e + "'")
		sys.exit(-1)

	# Launch calculation with current arguments
	print("Launching calculations with parameters '" + cur_args + "'")
	main(cur_args)
	



























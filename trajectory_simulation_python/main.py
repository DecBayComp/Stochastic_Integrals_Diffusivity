

## Dependencies
import argparse			# for command-line arguments
import csv 				# for csv files
import numpy as np
import os    			# for file operations
import random
import time 			# to measure elapsed time

from D_func import D_func
from constants import version, max_D_case, N as N_def, progress_update_interval, output_folder, str_mode, t_step, \
	internal_steps_number, CSV_DELIMITER, L, x_max, x_min


def main(arg_str):
	"""
	This is the new main file for simulations with ksi = alpha/D' as the only parameter.
	Use D_case = 2.
	This file replaces 'simulate_one_trajectory.py'
	Created: 2018-20-26
	"""
	
	# Random seed
	np.random.seed()


	## Define arguments
	arg_parser = argparse.ArgumentParser(description = 'Generation of a random walk trajectory')
	arg_parser.add_argument('-v', '--version', action = 'version', version = '%(prog)s ' + str(version))
	arg_parser.add_argument('-D', '--D_case', required = True, action = 'store', type = int, 
		choices = range(1, max_D_case + 1), metavar = 'D_case', help = 'D case to be simulated')
	arg_parser.add_argument('--id', required = True, action = 'store', type = int, 
		help = 'A simulation identifier to be added to the output file')
	arg_parser.add_argument('-N', action = 'store', type = float, help = 'Number of jumps to simulate')
	arg_parser.add_argument('--ksi', '--alpha_over_D', required = True, action = 'store', type = float, 
		help = 'Total force over the diffusivity gradient')


	## Analyze arguments
	input_args = arg_parser.parse_args(arg_str.split())


	## Use the analyzed arguments
	D_case = input_args.D_case
	file_id = input_args.id
	ksi = input_args.ksi

	# f_case = input_args.f_case
	# bl_random_lambda = input_args.rand
	# if bl_random_lambda:
		# Lambda = np.random.rand()
	# else:
		# Lambda = input_args.Lambda

	if input_args.N is not None:
		N = long(input_args.N)
	else:
		N=N_def
	update_progress_every = N / progress_update_interval


	# Confirm the parameters
	print('Calculating D_case: %i, alpha/D'': %.2f with N = %i steps' % (D_case, ksi, N))


	## Create output folder if missing
	if not os.path.isdir(output_folder):
		print('Output folder is missing. Creating')
		os.makedirs(output_folder)


	## Boundary conditions
	if str_mode == "periodic":
		bl_periodic = True
	else:
		bl_periodic = False
		

	## Initialize
	start_time = time.time()
	t_step_internal = t_step / internal_steps_number
	N_internal = N * internal_steps_number
	x_array = np.zeros(N+1, dtype = np.float32)
	dx_array = np.zeros(N, dtype = np.float32)
	t_mesh = np.arange(N + 1) * t_step


	## Perform a test save to ensure the process won't be killed
	print("Performing a test save to ensure enough memory is available...\n")
	filename = "sim_data_%09i.csv" % (file_id)
	output_full_path = output_folder + filename
	# print x_array
	output_data = np.zeros((N + 1, 2), dtype = np.float32)
	# print output_data
	# output_data = np.transpose(output_data)
	# Open the output file for writing
	with open(output_full_path, 'w') as file_pointer:
		csv_writer = csv.writer(file_pointer, delimiter = CSV_DELIMITER)
		csv_writer.writerows(output_data.tolist())
	# Remove the temporary file
	if os.path.isfile(output_full_path) :
		os.remove(output_full_path)
	print('Test file save succeeded. Proceeding to calculations\n')


	# Choosing the first point randomly
	x_0 = (np.random.rand() - 0.5) * L
	x_array[0] = x_0
	x_i = x_0
	# print x_i


	## Using Verlet method for iterations
	for i in range(N):
		q = np.random.normal(0.0, 1.0, internal_steps_number)    # Random noise
		dx_next = 0.0
		for m in range(internal_steps_number):
			# Calculate f and D at x_i
			[D_i, b_prime_b_i] = D_func(D_case, x_i, L)	# [D] = um^2/s, [D'] = [um/s]
			alpha_i = ksi * b_prime_b_i
					
			# Convert lists to np
			D_i = np.asarray(D_i)
			# b_prime_b_i = np.asarray(b_prime_b_i)
			b_i = np.sqrt(2.0 * D_i)
			
			alpha_i = np.asarray(alpha_i)
			
			# Create noise increment W_n (white noise)
			dW = np.sqrt(t_step_internal) * q[m]
			
			# Calculate increment (keping the terms up to dt)
			dx = alpha_i * t_step_internal + b_i * dW
					
			x_next = x_i + dx
			dx_next = dx_next + dx
			## Taking into account the BCs
			if bl_periodic:
				if x_next > x_max:
					x_next = x_next - L
				elif x_next < x_min:
					x_next = x_next + L
			else: 
				if x_next > x_max:
					x_next = 2.0 * x_max - x_next
				elif x_next < x_min:
					x_next = 2.0 * x_min - x_next
			# Save the starting point for the next internal round
			x_i = x_next      	        
		## Save
		x_array[i+1] = x_next
		dx_array[i] = dx_next
		# print x_i, x_next
		# print x_next
		
			
			 
		
		## Print out simulation progress
		if i % update_progress_every == 0:
			print("Simulation progress: %.1f %%. Elapsed time: %.2f s\n" % (float(i)/N*100.0, time.time() - start_time))


		# Iterating

	## Save data
	print("Saving trajectory...\n")

	# Prepare the output array
	output_data = np.zeros((N + 1, 2), dtype = np.float32)
	output_data[0, :] = [D_case, ksi]
	output_data[1:, 0] = x_array[0:N]
	output_data[1:, 1] = dx_array[0:N]


	# Open the output file for writing
	with open(output_full_path, 'w') as file_pointer:
		csv_writer = csv.writer(file_pointer, delimiter = CSV_DELIMITER, lineterminator = '\n')
		csv_writer.writerows(output_data.tolist())

	print('Trajectory saved successfully! Calculations finished in %.2f s\n' % (time.time() - start_time))
	print('Note: The first line in the ouput is the Lambda value.\n')









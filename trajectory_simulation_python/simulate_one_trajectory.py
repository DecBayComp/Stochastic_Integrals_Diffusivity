


## Dependencies
import argparse	# for command-line arguments
import csv 		# for csv files
import numpy as np
import os    # for file operations
import random
import time    # to measure elapsed time


## Constants
execfile("./constants.py")

## External functions
execfile("./D_func.py")
execfile("./f_func.py")

# Random seed
np.random.seed()


## Define arguments
arg_parser = argparse.ArgumentParser(description = 'Generation of a random walk trajectory')
arg_parser.add_argument('-v', '--version', action = 'version', version = '%(prog)s ' + str(version))
arg_parser.add_argument('-D', '--D_case', required = True, action = 'store', type = int, choices = range(1, max_D_case + 1), metavar = 'D_case', help = 'D case to be simulated')
arg_parser.add_argument('-f', '--f_case', required = True, action = 'store', type = int, choices = range(1, max_f_case + 1), metavar = 'f_case', help = 'f case to be simulated')
arg_parser.add_argument('--id', required = True, action = 'store', type = long, help = 'A simulation identifier to be added to the output file')
arg_group = arg_parser.add_mutually_exclusive_group(required = True)
arg_group.add_argument('-l', '--Lambda', type = float, help = 'The lambda* for simulations')
arg_group.add_argument('-r', '--rand', action = 'store_true', help = 'Use random lambda* for simulations')
arg_parser.add_argument('-N', action = 'store', type = float, help = 'Number of jumps to simulate')

## Analyze arguments
input_args = arg_parser.parse_args()

## Use the analyzed arguments
D_case = input_args.D_case
f_case = input_args.f_case
file_id = input_args.id
bl_random_lambda = input_args.rand
if bl_random_lambda:
	Lambda = np.random.rand();
	# output_folder = output_folder_rand 
else:
	Lambda = input_args.Lambda
	# output_folder = output_folder_fixed

if input_args.N is not None:
	N = long(input_args.N)
update_progress_every = N / 1.0e2;

# Confirm the parameters
print('Calculating D_case: %i, f_case: %i with N = %i steps' % (D_case, f_case, N))
# print input_args

## Create folder if missing
if not os.path.isdir(output_folder):
	print('Output folder is missing. Creating')
	os.makedirs(output_folder)


## Boundary conditions
if str_mode == "periodic":
	bl_periodic = True
else:
	bl_periodic = False

## Generate the output filename
# filename = 
	

## Initialize
start_time = time.time()
t_step_internal = t_step / internal_steps_number
N_internal = N * internal_steps_number
x_array = np.zeros(N+1, dtype = float)
# dx_array = np.zeros(N, dtype = float)
t_mesh = np.arange(N + 1) * t_step


## Perform a test save to ensure the process won't be killed
print("Performing a test save to ensure enough memory is available...\n")
filename = "sim_data_%09i.csv" % (file_id)
output_full_path = output_folder + filename
# print x_array
output_data = np.reshape(np.insert(x_array, 0, -1), (-1, 1))
# print output_data
# output_data = np.transpose(output_data)
# Open the output file for writing
with open(output_full_path, 'wb') as file_pointer:
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
	for m in range(internal_steps_number):
		# Calculate f and D at x_i
		[f_i, _] = f_func(f_case, x_i, L)
		[D_i, b_prime_b_i] = D_func(D_case, x_i, L)
		# Convert lists to arrays
		f_i = np.asarray(f_i)
		D_i = np.asarray(D_i)
		b_prime_b_i = np.asarray(b_prime_b_i)
		a_lambda_i = f_i * D_i / kBT
		b_i = np.sqrt(2.0 * D_i)
		
		# Create noise increment W_n
		dW = np.sqrt(t_step_internal) * q[m]
		
		# Calculate increment
		dx = (a_lambda_i * t_step_internal + b_i * dW +
			0.5 * b_prime_b_i * (dW**2 + (2.0 * Lambda - 1.0) * t_step_internal))
		# print dx
		
		x_next = x_i + dx
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
	# print x_i, x_next
	# print x_next
	
		
		 
	
	## Print out simulation progress
	if i % update_progress_every == 0:
		print("D case: %i. f case: %i. Simulation progress: %.1f %%. Elapsed time: %.2f s\n" % (D_case, f_case, float(i)/N*100.0, time.time() - start_time))
	# Iterating

## Save data
print("Saving trajectory...\n")
output_data = np.reshape(np.insert(x_array, 0, Lambda), (-1, 1))

# Open the output file for writing
with open(output_full_path, 'wb') as file_pointer:
	csv_writer = csv.writer(file_pointer, delimiter = CSV_DELIMITER)
	csv_writer.writerows(output_data.tolist())

print('Trajectory saved successfully! Calculations finished in %.2f s\n' % (time.time() - start_time))








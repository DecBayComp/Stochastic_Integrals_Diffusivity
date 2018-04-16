


## Constants
version = 20180411

dim = 2	# number of dimensions (1 or 2)
L = 1.0	# linear size of the system, in um 
x_min = 0.0	# in um
x_max = L	# in um
t_step = 0.04 # in seconds
N = int(1.0e5) # int(1.0e6) or int(1.0e5)


# Diffusivity constants
D_case = 1
D_0 = 0.01	# in um^2/s
D_ratio = 2.0	# Ratio between max and min of the diffusivity
q = 1	# number of D peaks
gamma_drag = 400.0	# viscous drag, in fN * s / um
k = 2.0 * q / L * (D_ratio - 1.0)	# D'/D
D_grad_abs = k * D_0	# absolute value of the diffusivity gradient used for force calculations


# Total force constants
ksi_range = (-10.0, 40.0)	# range of ratios of the total force to diffusivity gradient
ksi_step = 1.0	# sampling step in the interval


# Simulation parameters
progress_update_interval = 100.0
internal_steps_number = 100 # Integer. How many intermediate smaller steps are made before the next point is saved


# Batch parameters
trials = 1	# 1000
sleep_time = 0.2
logs_folder = "./logs/"
output_folder = "./output/"
args_file = "./arguments.dat"
args_lock = "./arguments.lock"
lock_timeout = 300

jobs_count_tars = 51	# 132+12
jobs_count_t_bayes = 132
jobs_count_onsager = 1
manager_script = "job_manager.py"
DETACHED_PROCESS = 0x00000008

output_folder = './output/';
CSV_DELIMITER = ';';


# k = 100.0 # diffusivity gradient, D'/D, in um^-1

max_D_case = 7;

max_f_case = 8;

str_mode = 'periodic'; #    bc_type = ENUM_BC_PERIODIC;













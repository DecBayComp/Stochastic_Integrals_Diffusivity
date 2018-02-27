


## Constants
version = 20180226


# % m = 1;
# kBT = 1.0
L = 1.0	# in um
x_min = -L/2.0	# in um
x_max = L/2.0	# in um
gamma_drag = 400.0	# viscous drag, in fN * s / um
# % % L = x_max - x_min;
# % bl_periodic = true;
# % T = 100000;
t_step = 0.125 # in seconds
N = int(1.0e5) # int(1.0e6)
progress_update_interval = 100.0
internal_steps_number = 100 # Integer. How many intermediate smaller steps are made before the next point is saved


## Batch parameters
trials = 20
sleep_time = 0.2
logs_folder = "./logs/"
output_folder = "./output/"
args_file = "./arguments.dat"
ksi_range = (-1, 2)
ksi_step = 0.1


output_folder = './output/';
CSV_DELIMITER = ';';



## Diffusivity constants
D_0 = 0.01	# in um^2/s
k = 2.0 # in um^-1
max_D_case = 7;



max_f_case = 8;


# Choosing the boundary conditions
# str_mode = "inf_walls" #     bc_type = ENUM_BC_INF_WALLS;
str_mode = 'periodic'; #    bc_type = ENUM_BC_PERIODIC;













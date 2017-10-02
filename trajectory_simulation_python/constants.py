


## Constants
version = 1.01
lambda_array = [0.0, 0.5, 1.0]
lambda_count = len(lambda_array)


# % m = 1;
kBT = 1.0
L = 1.0	# in um
x_min = -L/2.0	# in um
x_max = L/2.0	# in um
gamma_drag = 400.0	# viscous drag, in fN * s / um
# % % L = x_max - x_min;
# % bl_periodic = true;
# % T = 100000;
t_step = 1.0e-2 # in seconds, L^2 / D_max / 100 = 10 / 1 / 100 = 0.1
N = int(1e5) 
internal_steps_number = 100 # Integer. How many intermediate smaller steps are made before the next point is saved
# SIGMA_MIN = 0;
# SIMGA_MAX = 10;
# D_MIN = 0;
# D_MAX = 100;
# ABS_TOLERANCE = 1e-7;
# REL_TOLERANCE = 1e-7;
# % selected_x_over_L_coordinates = [-0.25, 0, 0.25];
selected_x_over_L = 0.25
# fD_marginalized_steps = 1 + 2^3;    % 1 + 2^5
# min_points_in_bin = 10;
# SIMULATION_TRIES_PER_CASE = 30;
# bl_use_adaptive_mesh = true;



# marker_size = 7;
# font_size = 16;

# output_figures_folder = './figures_for_article/';
# output_data_folder = './processed_data/';
output_folder = './output/';
# output_folder_rand = './output/rand/';
# output_folder_fixed = './output/fixed/';
CSV_DELIMITER = ';';
# bl_save_figures = true;
# bl_save_data = true;
# % bl_save_data = false;


max_D_case = 6;
max_f_case = 7;


# %% Choosing the boundary conditions
# ENUM_BC_INF_WALLS = 1;
# ENUM_BC_PERIODIC = 2;
str_mode = "inf_walls" #     bc_type = ENUM_BC_INF_WALLS;
# str_mode = 'periodic'; #    bc_type = ENUM_BC_PERIODIC;



# lambda_names_array = {'Ito', 'Stratonovich', 'Isothermal'}; 







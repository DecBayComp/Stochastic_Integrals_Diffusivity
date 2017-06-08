
%% Constants
% m = 1;
kBT = 1;
L = 1;
x_min = -L/2;
x_max = L/2;
% % L = x_max - x_min;
% bl_periodic = true;
% T = 100000;
t_step = 1e-2;   % L^2 / D_max / 100 = 10 / 1 / 100 = 0.1
N = 1e5;        % Times the system was explored: (N*t_step) / (L^2 / D_max) = N*t_step*D_max/L^2 = 1e4*0.1*1/100 = 10   OK
min_points_in_bin = 500;
T = t_step * N;
internal_steps_number = 100;    % Integer. How many intermediate smaller steps are made before the next point is saved
update_progress_every = 100;
SIGMA_MIN = 0;
SIMGA_MAX = 10;
D_MIN = 0;
D_MAX = 100;
ABS_TOLERANCE = 1e-7;
REL_TOLERANCE = 1e-7;
% selected_x_over_L_coordinates = [-0.25, 0, 0.25];
selected_x_over_L = 0.4;    % 0.4
fD_marginalized_steps = 1 + 2^5;    % 1 + 2^5
SIMULATION_TRIES_PER_CASE = 30;
bl_use_adaptive_mesh = true;
alpha_smoooth = 1e-3;
fine_mesh_steps_count = 1000 + 1;
CONF_LEVEL = 0.95;


%% Regularization parameters
alpha_reg = 1e0;


%% Plot parameters
marker_size = 6;
font_size = 12;
legend_font_size = 9;
markers_list = {'x', 'o', '+', '^', 's','d','v'};
line_width = 1.2;
% markers_list = {'-o','-s','-d','-^','-v'};



output_figures_folder = './figures_for_article/';
output_data_folder = './processed_data/';
input_data_folder = 'D:\Calculations_data\ito-stratonovich\output\';
CSV_DELIMITER = ';';
bl_save_figures = true;
bl_save_data = true;
% bl_save_data = false;


max_D_case_number = 6;
max_f_case_number = 7;


%% Plotting results for the article
selected_D_case = 6;
selected_f_case = 7;


%% Choosing the boundary conditions
ENUM_BC_INF_WALLS = 1;
ENUM_BC_PERIODIC = 2;
str_mode = 'inf_walls';     bc_type = ENUM_BC_INF_WALLS;
% str_mode = 'periodic';    bc_type = ENUM_BC_PERIODIC;


% lambda_array = [0, 0.5, 1];
lambda_names_array = {'Ito�', 'Stratonovich', 'Isothermal'}; 
% lambda_count = length(lambda_array);


% if N>1e4
%     bl_short_trajectory = false;
% else
    bl_short_trajectory = true;
% end;


%% Enumerate conventions
enum_conv_Ito = 1;
enum_conv_Stratonovich = 2;
enum_conv_Hanggi = 3;
enum_conv_divine = 4;
enum_conv_marginalized = 5;
conventions_count = 5;
conventions_names = {'Ito', 'Str', 'Hng', 'Orcl', 'Mar'};


%% Enumerate lambda simulation types
enum_lambda_Ito = 1;
enum_lambda_Stratonovich = 2;
enum_lambda_Hanggi = 3;
enum_lambda_rand = 4;
lambda_types_count = 4;
lambda_types_names = {'Ito', 'Str', 'Hng', 'rnd'};


% Define colors
color_sequence = [0    0.4470    0.7410;...
                    0.9290    0.6940    0.1250;...
                    0.8500    0.3250    0.0980;...
                    139/255 87/255 66/255;...   %	Brown               
                    107/255 142/255 35/255;...   %	Olive                    
                    0.502 0.502 1;...   %	Light purple                    
                    0 0.251 0;...   %	Dark green       
                    0 0 1;...   %	Blue
                    0.9412 0.4706 0;... %   Orange
                    0.502 0.251 0;...   %	Brown                    
                    0 0.502 0.502;...   %	Turquoise
                    1 0 0;...   %	Bright red
%                     1 1 1;...   %	White
                    1 0.502 0.502;...   %	Peach
                    0 1 1;...   %	Cyan
%                     0.502 0.502 0.502;...   %	Gray
                    0.502 0 0;...   %	Burgundy 
                    1 0 1;...   %	Pink
                    0.251 0 0.502;...   %	Purple
                    1 1 0;...  %	Yellow                    
                    0 0 0;...   %	Black
                    0 0.251 0;...   %	Dark green       
                    0 0 1;...   %	Blue
                    0.9412 0.4706 0;... %   Orange
                    0.502 0.251 0;...   %	Brown                    
                    0 0.502 0.502;...   %	Turquoise
                    1 0 0;...   %	Bright red
%                     1 1 1;...   %	White
                    1 0.502 0.502;...   %	Peach
                    0 1 1;...   %	Cyan
%                     0.502 0.502 0.502;...   %	Gray
                    0 1 0;...   %	Bright green
                    1 0 1;...   %	Pink
                    0.251 0 0.502;...   %	Purple                    
                    1 1 0;...  %	Yellow                    
                    0 0 0];...   %	Black
                    
my_constructs_color_sequence = [ 0.502 0 0;...   %	Burgundy 
                            1 0 1;...   %	Pink
                            0.251 0 0.502;...   %	Purple
                            1 1 0];  %	Yellow



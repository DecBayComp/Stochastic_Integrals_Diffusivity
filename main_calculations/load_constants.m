
%% Constants
% m = 1;
kBT = 1;	% change to 4.14 * 10^(-21) J
L = 1;	% in um
x_min = -L/2;
x_max = L/2;
gamma_drag = 400.0;	% viscous drag, in fN * s / um
% % L = x_max - x_min;
% bl_periodic = true;
% T = 100000;
t_step = 0.08;   % L^2 / D_max / 100 = 10 / 1 / 100 = 0.1
N = 5.0e5;        % Times the system was explored: (N*t_step) / (L^2 / D_max) = N*t_step*D_max/L^2 = 1e4*0.1*1/100 = 10   OK
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
fD_marginalized_steps = 1 + 2^5;    % 1 + 2^5
SIMULATION_TRIES_PER_CASE = 30;
bl_use_adaptive_mesh = true;
alpha_smoooth = 1e-3;
fine_mesh_steps_count = 1000 + 1;
CONF_LEVEL = 0.95;
w = 10;


%% Binning
points_in_bin_avg = 1e4;
min_bin_to_jump_ratio = 2;	% require a bin to be at least several times larger than the mean jump in it. This corresponds to a 61% probability to stay in bin after jump
bl_keep_only_min_points_in_bin = true;	% When true, only the minimum number of points is kept per bin. Extra points are randomly omitted within a bin


%% Regularization parameters
alpha_reg = 0.1;


%% Plot parameters
marker_size = 6;
font_size = 12;
subplot_label_font_size = 16;
legend_font_size = font_size - 3;
markers_list = {'s', 'o', '^', '+', 'x', 'd','v'};
line_width = 1.2;
line_width_theor = line_width - 0.5;

% Bin colors
load_color_scheme;
bin_color = my_colors(3).White;

% markers_list = {'-o','-s','-d','-^','-v'};



output_figures_folder = './figures_for_article/';
output_data_folder = './processed_data/';
input_data_folder = '/home/aserov/Documents/Calculated_data/two_forces/';  % Ubuntu
% input_data_folder = '/home/aserov/Documents/Calculated_data/dilemma_no_force/';  % Ubuntu
% input_data_folder = '/Users/alexander_serov/Calculations_data/ito-stratonovich/'; % Mac
fail_rates_filename = 'Fail_rates.dat';
CSV_DELIMITER = ';';
bl_save_figures = true;
bl_save_data = true;
% bl_save_data = false;


max_D_case_number = 6;
max_f_case_number = 8;


%% Plotting results for the article
selected_D_case = 2;
selected_f_case = 8;


%% Choosing the boundary conditions
ENUM_BC_INF_WALLS = 1;
ENUM_BC_PERIODIC = 2;
str_mode = 'inf_walls';     bc_type = ENUM_BC_INF_WALLS;
% str_mode = 'periodic';    bc_type = ENUM_BC_PERIODIC;


% lambda_array = [0, 0.5, 1];
lambda_names_array = {'Ito', 'Stratonovich', 'Isothermal'}; 
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
conventions_tex_names = {'It\^o', 'Stratonovich', 'H\"anggi', 'Oracle', 'Marginalized'};


%% Enumerate lambda simulation types
enum_lambda_Ito = 1;
enum_lambda_Stratonovich = 2;
enum_lambda_Hanggi = 3;
enum_lambda_rand = 4;
lambda_types_count = 4;
lambda_types_names = {'Ito', 'Str', 'Hng', 'rnd'};
lambda_types_tex_names = {'It\^o', 'Stratonovich', 'H\"anggi', 'Random'};
lambda_ind_for_KS_plot = 3;


% Define colors
load_color_scheme;
color_sequence = [standard_colors(1).DeepBlue; my_colors(5).Green; my_colors(1).Orange; my_colors(1).WarmBrown; standard_colors(1).Purple];
% color_sequence = [0    0.4470    0.7410;...
%                     0.9290    0.6940    0.1250;...
%                     0.8500    0.3250    0.0980;...
%                     139/255 87/255 66/255;...   %	Brown               
%                     107/255 142/255 35/255;...   %	Olive                    
%                     0.502 0.502 1;...   %	Light purple                    
%                     0 0.251 0;...   %	Dark green       
%                     0 0 1;...   %	Blue
%                     0.9412 0.4706 0;... %   Orange
%                     0.502 0.251 0;...   %	Brown                    
%                     0 0.502 0.502;...   %	Turquoise
%                     1 0 0;...   %	Bright red
% %                     1 1 1;...   %	White
%                     1 0.502 0.502;...   %	Peach
%                     0 1 1;...   %	Cyan
% %                     0.502 0.502 0.502;...   %	Gray
%                     0.502 0 0;...   %	Burgundy 
%                     1 0 1;...   %	Pink
%                     0.251 0 0.502;...   %	Purple
%                     1 1 0;...  %	Yellow                    
%                     0 0 0;...   %	Black
%                     0 0.251 0;...   %	Dark green       
%                     0 0 1;...   %	Blue
%                     0.9412 0.4706 0;... %   Orange
%                     0.502 0.251 0;...   %	Brown                    
%                     0 0.502 0.502;...   %	Turquoise
%                     1 0 0;...   %	Bright red
% %                     1 1 1;...   %	White
%                     1 0.502 0.502;...   %	Peach
%                     0 1 1;...   %	Cyan
% %                     0.502 0.502 0.502;...   %	Gray
%                     0 1 0;...   %	Bright green
%                     1 0 1;...   %	Pink
%                     0.251 0 0.502;...   %	Purple                    
%                     1 1 0;...  %	Yellow                    
%                     0 0 0];...   %	Black
                    
my_constructs_color_sequence = [ 0.502 0 0;...   %	Burgundy 
                            1 0 1;...   %	Pink
                            0.251 0 0.502;...   %	Purple
                            1 1 0];  %	Yellow




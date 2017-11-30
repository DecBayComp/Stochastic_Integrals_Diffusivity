%% Plot histograms of the global Bayes factor for different simulation and inference conventions



function plot_article_global_Bayes_factor(data_struct, trials_data, fig_count, bl_save_figures)



%% Constants
load_constants;

% Subplot parameters
SH = 0.05;
SV = 0.12;
ML = 0.07;
MR = 0.02;
MT = 0.07;
MB = 0.17;
rows = 1;
cols = 4;

face_alpha = 0.5;



%% Initialize
% Initalize plot
h_fig = figure(fig_count);
set_article_figure_size(h_fig, rows, 2, 1);
clf;

% Initalize subplots
subaxis(rows, cols, 1, 'SH', SH, 'SV', SV, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);

% Load Bayes factor data. Format: [trial x convention]
trial_simulation_type = data_struct.trial_simulation_type;
log_K_G_all_trials = data_struct.trials_log_K_G;


% Plot
subaxis(1);
hold on;
lambda_type = enum_lambda_Ito;

for convention = 1:conventions_count
	% Extract current lambda type and convention data for log_K_G
	log_K_G = log_K_G_all_trials(trial_simulation_type == lambda_type, convention);

	% Plot
	histogram(log_K_G, 'Normalization', 'probability', 'FaceColor', color_sequence(convention, :));
	
end

xlabel('$\ln K_G$', 'interpreter', 'latex');





1;






















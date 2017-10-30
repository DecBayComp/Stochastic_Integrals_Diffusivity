%% Plot MAP distribution of b and posterior in one bin. Plot the cumulative distributions and the KS distance between the two


function plot_article_b_KS_distance(data_struct, trials_data, fig_count, bl_save_figures)


%% Constants
load_constants;
% b_min = 0.1;
% b_max = 0.2;
b_steps = 1e2;
b_extend_factor = 1.5;
REL_TOLERANCE = 1e-3;
ABS_TOLERANCE = 1e-3;

trials_to_plot = 5;	% Number of trials' priors to plot

height_factor = 1;
output_filename = 'b_MAP_vs_posteriors.pdf';

% % Subplot params
% rows = 2;
% cols = 1;
% SH = 0.08;
% SV = 0.1;
% ML = 0.12;
% MR = 0.0125;
% MT = 0.08;
% MB = 0.15;




%% Calculate the empirical distributions
% Identify the bin with the selected x coordinate (same for all trials and lambda)
[selected_bins_indices, selected_bins_centers] = get_selected_bins_indices(data_struct);
bin = selected_bins_indices;

% Load MAP distribution for current lambda type
trial_simulation_type = data_struct.trial_simulation_type;
b_MAP_data = data_struct.trials_MAP_b(trial_simulation_type == lambda_ind_for_KS_plot, bin, 1);

% Drop NaN
b_MAP_data = b_MAP_data(~isnan(b_MAP_data));
n = length(b_MAP_data);

% Sort
b_MAP_data = sort(b_MAP_data);

% Calculate b_min and b_max
b_min = min(b_MAP_data);
b_max = max(b_MAP_data);

% Extend the interval by 2 times
b_width = b_max - b_min;
b_max = b_max + b_width/2 * (b_extend_factor - 1);
b_min = b_min - b_width/2 * (b_extend_factor - 1);

% Calculate b_step
b_step = (b_max - b_min) / (b_steps - 1);



%% Calculate empirical cdf for the MAP values
% Initialize
n_double_mesh = zeros(2 * n, 1);
b_MAP_cdf = zeros(2 * n, 1);

% n_double_mesh(1:2:2 * n-1) = 1:n;
% n_double_mesh(2:2:2 * n) = 1:n;

% Double b mesh for plot. The double mesh is need to plot discontinuous jumps
b_double_mesh = zeros(2 * n, 1);
b_double_mesh(1:2:2 * n-1) = b_MAP_data;
b_double_mesh(2:2:2 * n) = b_MAP_data;

% Calculate CDF on the double b mesh
b_MAP_cdf(1:2:2 * n-1) = (0:n-1)/n;
b_MAP_cdf(2:2:2 * n) = (1:n)/n;

% % Check the KS distance
% calculate_KS_distance(b_MAP_data, b_posterior_func_wrap)



%% Get b posterior in the selected bin
% Initialize
b_posterior_pdf = zeros(trials_to_plot, b_steps);
b_posterior_cdf = zeros(trials_to_plot, b_steps);

% Trials to consider
first_trial = data_struct.trial_first_simulation_type_index(lambda_ind_for_KS_plot);
trials_indices = first_trial + (0:trials_to_plot-1) * lambda_types_count;

% Load data
for k = 1:trials_to_plot
	trial = trials_indices(k);
	cur_data_struct = trials_data{trial};

	% Load the posterior for this bin and trial
	b_posterior_func_wrap = @(b) bin_b_posterior_func (bin, b, t_step, cur_data_struct, 'forward');

	% Define the continuous CDF
	b_posterior_cdf_func_wrap = @(b) integral(b_posterior_func_wrap, -Inf, b, 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE);

	% Prepare b mesh
	b_mesh = b_min:b_step:b_max;

	% Calculate b posterior pdf
	b_posterior_pdf(k, :) = b_posterior_func_wrap(b_mesh);

	% Calculate continuous cdf
	for i = 1:b_steps
		b_posterior_cdf(k, i) = b_posterior_cdf_func_wrap(b_mesh(i));
	end;
	
end;



%% Initialize figure
h_fig = figure(fig_count);
set_article_figure_size(h_fig, 1, 1, height_factor);
clf;



%% == (A): PDF ==
% Initialize subplot
% h_sub = subaxis(rows, cols, 1, 'SH', SH, 'SV', SV, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);
hold on;

% Plot histogram of MAP b
histogram(b_MAP_data, 'Normalization', 'pdf', 'FaceColor', color_sequence(2, :));

%% Plot b posterior
plot(b_mesh, b_posterior_pdf, 'LineWidth', line_width, 'color', color_sequence(1, :));


%% Adjust
xlabel('$b$, $\mu \mathrm{m} / s^{1/2}$', 'interpreter', 'latex');
ylabel('PDF', 'interpreter', 'latex');
tmp_str = sprintf('$\\lambda^* = 1.0$, $x \\approx %.2f~\\mu \\mathrm{m}$', selected_bins_centers);
title(tmp_str, 'interpreter', 'latex');
grid on;
box on;



% % % %% == (B): CDF ==
% % % % Initialize subplot
% % % h_sub = subaxis(rows, cols, 2);
% % % hold on;
% % % 
% % % % Plot the bin posterior
% % % plot(b_mesh, b_posterior_cdf, 'r', 'LineWidth', line_width, 'color', color_sequence(1, :));
% % % 
% % % % Plot the empirical CDF
% % % plot_step = 100;
% % % plot(b_double_mesh(1:plot_step:end), b_MAP_cdf(1:plot_step:end), '-x', 'LineWidth', line_width, 'color', color_sequence(2, :));
% % % 
% % % 
% % % %% Adjust
% % % xlabel('$b$, $\mu m / s^{1/2}$', 'interpreter', 'latex');
% % % 
% % % % tmp_str = sprintf('$\\lambda^* = 1.0$, $x \\approx %.2f~\\mu m$', selected_bins_centers);
% % % % title(tmp_str, 'interpreter', 'latex');
% % % grid on;
% % % ylim([-0.01,1.01]);



%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize', [fig_pos(3), fig_pos(4)]);

% Set filename
output_full_path = strcat(output_figures_folder, output_filename);
if bl_save_figures
    print(h_fig, output_full_path, '-dpdf', '-r0');
end;









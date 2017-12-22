%% Plot Bayes factor profile for each simulation and inference model combination


function plot_article_local_bayes_factor(data_struct, bl_force, fig_count, bl_save_figures)



%% Constants
load_constants;
sublabel_x = 0.03;
sublabel_y = 0.08;
x_lim_vec = [x_min, x_max];
y_lim_vec = [-1, 1] * 1.3;
y_lim_vec_profile = [-1, 1] * 0.18;
output_filename_base = 'K_L';
% Subplot parameters
SH = 0.03;
SV = 0.12;
ML = 0.06;
MR = 0.02;
MT = 0.1;
MB = 0.15;
rows = 1;
cols = 4;
jitter_scale = 1/6;

% Skip some bins
plot_every = 1;

% Label params
sublabel_x = 0.015;
sublabel_y = 1.12;

x_tick_increment = 0.2;

% Confidence zones
% CONF_LN_K_POSITIVE = [1; 3];
% CONF_LN_K_STRONG = [3; 5];
% CONF_LN_K_VERY_STRONG = [5; 100];
CONF_LN_K_STRONG = 3;

% Set confidence zones colors
i=3;
conf_color_sequence = [my_colors(i).White; my_colors(i+1).White; my_colors(i+2).White];
% conf_color_sequence = [233, 243, 216; 234, 246, 255]/255;


%% Load data
% log_K_L_mean = data_struct.log_K_L_mean;
x_bins_centers = data_struct.x_bins_centers;
bins_count = length(x_bins_centers);
trial_simulation_type = data_struct.trial_simulation_type;
bin_widths = data_struct.x_bins_widths;
bin_borders = [1; 1] * data_struct.x_bins_centers' + [-1/2; 1/2] * bin_widths';
log_K_L_all_trials = data_struct.trials_log_K_L;	% [trials x bins x conventions]

% Generate normally distributed random jitter
% jitter = bin_widths * randn(1, conventions_count) * jitter_scale;
jitter = ones(bins_count, 1) * ((0:conventions_count-1) / (conventions_count - 1) - 1/2) * 2 * jitter_scale;



%% Plot
% Initalize plot
h_fig = figure(fig_count);
set_article_figure_size(h_fig, rows, 2, 1);
clf;

% Initalize subplots
subaxis(rows, cols, 1, 'SH', SH, 'SV', SV, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);
% for lambda_type = 1:lambda_types_count



%% Calculate
mean_log_K_L = zeros(lambda_types_count, bins_count, conventions_count);
std_log_K_L = zeros(lambda_types_count, bins_count, conventions_count);
eb_log_K_L = zeros(lambda_types_count, bins_count, conventions_count);	% error bars
mean_evidence = zeros(lambda_types_count, bins_count, conventions_count);
for lambda_type = 1:lambda_types_count
	% Extract current-lambda-type results
	log_K_L = log_K_L_all_trials(trial_simulation_type == lambda_type, :, :);
	
	% Calculate mean and std of the log(K)! over all trials of the same type
	mean_log_K_L(lambda_type, :, :) = mean(log_K_L, 1);
	std_log_K_L(lambda_type, :, :) = std(log_K_L, [], 1);
	eb_log_K_L(lambda_type, :, :) = std_log_K_L(lambda_type, :, :) * sqrt(2) * erfinv(0.95);
end

% Evaluate evidence for the presence of force on 3-values scale
mean_evidence(mean_log_K_L >= CONF_LN_K_STRONG) = 1;
mean_evidence(mean_log_K_L <= -CONF_LN_K_STRONG) = -1;


	
for lambda_type = 1:lambda_types_count	
	%% Plot
	subaxis(lambda_type);
	hold on;
	str_legend = {};
	for convention = 1:conventions_count
		plot(x_bins_centers, squeeze(mean_evidence(lambda_type, :, convention)) + jitter(:, convention)', strcat('', markers_list{convention}),...
			'color', color_sequence(convention, :), 'markers', marker_size, 'linewidth', line_width);
		str_legend{length(str_legend) + 1} = conventions_names{convention};
	end;
% 	x_lim_vec = xlim();
	xlim(x_lim_vec);

	% Theory
	h_theor = plot(x_lim_vec, [0,0], 'k--', 'linewidth', line_width_theor);

	% Adjust
	ylim(y_lim_vec);
% 	ylim('auto');

	% Subplot label
	text(sublabel_x, sublabel_y, char('A' + lambda_type - 1 + 4 * bl_force), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', subplot_label_font_size);

	% Modify ticks
	set(gca,'xtick', x_min:x_tick_increment:x_max);
	set(gca,'ytick', [])
	if lambda_type == 1
		set(gca,'ytick', [-1, 0, 1]);
	end;


	% Axes labels
	xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
	if lambda_type == 1 && ~bl_force
% 		ylabel('Bayes factor $\langle \ln K_L \rangle$', 'interpreter', 'latex');
		ylabel('$\langle \ln K_L \rangle$ for spurious-force model', 'interpreter', 'latex');
	elseif lambda_type == 1 && bl_force
		ylabel('$\langle \ln K_L \rangle$ for local-force model', 'interpreter', 'latex');
	end

	% Title
	str_title = sprintf('%s sim.', lambda_types_tex_names{lambda_type});
	title(str_title, 'interpreter', 'latex');
	
	% Legend
	if lambda_type == 2 && ~data_struct.bl_force
		legend(str_legend, 'location', 'northeast', 'fontsize', legend_font_size);
	
	elseif lambda_type == 2 && data_struct.bl_force
		legend(str_legend, 'location', 'southwest', 'fontsize', legend_font_size);
	end;
	
	% Add confidence zones
	x_lim_vec = xlim();
	
	% Color bin borders
	color_bins(bin_borders, y_lim_vec, bin_color);
	
	% Put axes on top
	set(gca, 'Layer', 'top');

end;



%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize', [fig_pos(3), fig_pos(4)]);
% Set filename
output_filename = strcat(output_filename_base, '_', data_struct.str_force, '.pdf');
output_full_path = strcat(output_figures_folder, output_filename);
if bl_save_figures
    print(h_fig, output_full_path, '-dpdf', '-r0');
end;










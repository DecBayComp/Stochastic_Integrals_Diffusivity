%% Plot histograms of the global Bayes factor for different simulation and inference conventions



function plot_article_global_bayes_factor(data_struct, fig_count, bl_save_figures)



%% Constants
load_constants;

% Subplot parameters
SH = 0.05;
SV = 0.12;
ML = 0.07;
MR = 0.02;
MT = 0.09;
MB = 0.17;
rows = 1;
cols = 4;

% Label params
sublabel_x = 0.015;
sublabel_y = 1.12;

y_increase_frac = 0.05;



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

% Calculate std of the global Bayes factor for each simulation - inference couple
log_mean_K_G = zeros(lambda_types_count, conventions_count);
log_std_K_G = zeros(lambda_types_count, conventions_count);
y_lim_vec = [0,0];
for lambda_type = 1:lambda_types_count
	% Initialize plot
	subaxis(lambda_type);
	hold on;
	
	% Extract current lambda type results
	log_K_G = log_K_G_all_trials(trial_simulation_type == lambda_type, :);
	
	% Calculate mean
	log_mean_K_G(lambda_type, :) = log(mean(exp(log_K_G), 1));
	
	% Calculate stds of the global Bayes factor
% 	for convention = 1:conventions_count
	K_G_std(lambda_type, :) = std(exp(log_K_G), [], 1);
	log_std_K_G(lambda_type, :) = log(std(exp(log_K_G), [], 1));
% 	end
	
	% Plot with error bars
	str_legend = {};
	for convention = 1:conventions_count
		errorbar(convention, log_mean_K_G(lambda_type, convention), log_std_K_G(lambda_type, convention) * sqrt(2) * erfinv(0.95), markers_list{convention},...
			'color', color_sequence(convention, :), 'markers', marker_size, 'linewidth', line_width);
		str_legend{length(str_legend) + 1} = conventions_names{convention};
	end;
	
	% Add x labels
	set(gca,'XTick',1:conventions_count, 'XTickLabel', conventions_names);
	
	% Subplot label
	text(sublabel_x, sublabel_y, char('A' + lambda_type - 1), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', subplot_label_font_size);

	% Axes labels
	if lambda_type ==1
		ylabel('Global Bayes factor $\ln \langle K_G \rangle$', 'interpreter', 'latex');
	end;

	% Title
	str_title = sprintf('%s sim.', lambda_types_tex_names{lambda_type});
	title(str_title, 'interpreter', 'latex');
	
	% Theory
	h_theor = plot(xlim(), [0,0], 'k--', 'linewidth', line_width_theor);
	uistack(h_theor, 'bottom');
	
% % % 	% Keep y limits if larger
% % % 	cur_y_lim_vec = ylim();
% % % 	y_lim_vec(1) = min(y_lim_vec(1), cur_y_lim_vec(1));
% % % 	y_lim_vec(2) = max(y_lim_vec(2), cur_y_lim_vec(2));
	
	
% % % 	% Legend
% % % 	if lambda_type == 1
% % % 		legend(str_legend, 'location', 'northeast');
% % % 	end;
	
end



%% Reset y limits
y_lim_vec = [0, 0];
y_lim_vec(1) = min(min(log_mean_K_G - log_std_K_G * sqrt(2) * erfinv(0.95)));
y_lim_vec(2) = max(max(log_mean_K_G + log_std_K_G * sqrt(2) * erfinv(0.95)));

% Slightly increase the interval width
y_lim_vec = y_lim_vec + [-1/2, 1/2]  * (y_lim_vec(2) - y_lim_vec(1)) * y_increase_frac;


for lambda_type = 1:lambda_types_count
	subaxis(lambda_type);
	ylim(y_lim_vec);
end;

% % % 
% % % % Plot
% % % 
% % % lambda_type = enum_lambda_Ito;
% % % 
% % % for convention = 1:conventions_count
% % % 	% Extract current lambda type and convention data for log_K_G
% % % 	log_K_G = log_K_G_all_trials(trial_simulation_type == lambda_type, convention);
% % % 
% % % 	% Plot
% % % % 	histogram(log_K_G, 'Normalization', 'probability', 'FaceColor', color_sequence(convention, :));
% % % 	
% % % end
% % % 
% % % xlabel('$\ln K_G$', 'interpreter', 'latex');





1;






















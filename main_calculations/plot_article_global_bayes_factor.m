%% Plot histograms of the global Bayes factor for different simulation and inference conventions



function plot_article_global_bayes_factor(data_struct, bl_force, fig_count, bl_save_figures)



%% Constants
load_constants;
output_filename_base = 'K_G';

% Overrides
line_width = 1.25;
marker_size = marker_size + 1;

% Figure size parameters
page_width_frac = 1;
height_factor = 0.7;

% Subplot parameters
SH = 0.04;
SV = 0.12;
ML = 0.075;
MR = 0.005;
MT = 0.1;
MB = 0.12;
rows = 1;
cols = 4;

% Label params
sublabel_x = 0;
sublabel_y = 0.8 / height_factor;

y_increase_frac = 0.05;



%% Initialize
% Initalize plot
h_fig = figure(fig_count);
set_article_figure_size(h_fig, rows, page_width_frac, height_factor);
clf;

% Initalize subplots
subaxis(rows, cols, 1, 'SH', SH, 'SV', SV, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);

% % Set initial sublabel
% if ~bl_force
% 	chr_sublabel_start = 'A';
% else
% 	chr_sublabel_start = '';
% end

% Load Bayes factor data. Format: [trial x convention]
trial_simulation_type = data_struct.trial_simulation_type;
log_K_G_all_trials = data_struct.trials_log_K_G;

% Initialize
mean_log_K_G = zeros(lambda_types_count, conventions_count);
std_log_K_G = zeros(lambda_types_count, conventions_count);
eb_log_K_G = zeros(lambda_types_count, conventions_count);	% error bars
for lambda_type = 1:lambda_types_count
	% Initialize plot
	subaxis(lambda_type);
	hold on;
	
	% Extract current-lambda-type results [trials x conventions]
	log_K_G = log_K_G_all_trials(trial_simulation_type == lambda_type, :);
	
	% Calculate mean and std of the log(K)!
	mean_log_K_G(lambda_type, :) = mean(log_K_G, 1);
	std_log_K_G(lambda_type, :) = std(log_K_G, [], 1);
	eb_log_K_G(lambda_type, :) = std_log_K_G(lambda_type, :) * sqrt(2) * erfinv(0.95);
	
	% Plot with error bars
	str_legend = cell(conventions_count, 1);
	for convention = 1:conventions_count
		errorbar(convention, mean_log_K_G(lambda_type, convention), eb_log_K_G(lambda_type, convention), markers_list{convention},...
			'color', color_sequence(convention, :), 'markers', marker_size, 'linewidth', line_width);
		str_legend{convention} = conventions_names{convention};
	end
	
	% Add x labels
	set(gca,'XTick',1:conventions_count, 'XTickLabel', conventions_names, 'FontSize', font_size);
    xlim([0.5, conventions_count + 0.5]);
	
	% Subplot label
	text(sublabel_x, sublabel_y, char('A' + lambda_type - 1 + 4 * bl_force), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', subplot_label_font_size);

	% Axes labels
	if lambda_type == 1 && ~bl_force
% 		ylabel('Bayes factor $\langle \ln K_L \rangle$', 'interpreter', 'latex');
		ylabel({'No force', '$\langle \ln K_G \rangle$'}, 'interpreter', 'latex');
	elseif lambda_type == 1 && bl_force
		ylabel({'Force', '$\langle \ln K_G \rangle$'}, 'interpreter', 'latex');
    end
    
    

	% Title
	str_title = sprintf('%s sim.', lambda_types_tex_names{lambda_type});
	title(str_title, 'interpreter', 'latex');
	
	% Theory
	h_theor = plot(xlim(), [0,0], '-', 'linewidth', line_width_theor, 'color', axes_color);
	uistack(h_theor, 'bottom');
	
% % % 	% Keep y limits if larger
% % % 	cur_y_lim_vec = ylim();
% % % 	y_lim_vec(1) = min(y_lim_vec(1), cur_y_lim_vec(1));
% % % 	y_lim_vec(2) = max(y_lim_vec(2), cur_y_lim_vec(2));
	
	
% 	% Legend
% 	if lambda_type == 1
% 		legend(str_legend, 'location', 'northeast');
%     end
	
end



%% Reset y limits
y_lim_vec = [0, 0];
y_lim_vec(1) = min(min(mean_log_K_G - eb_log_K_G));
y_lim_vec(2) = max(max(mean_log_K_G + eb_log_K_G));

% Slightly increase the interval width
y_lim_vec = y_lim_vec + [-1/2, 1/2]  * (y_lim_vec(2) - y_lim_vec(1)) * y_increase_frac;

% Model-specific modifications
if bl_force
   y_lim_vec(1) = 0;
end


for lambda_type = 1:lambda_types_count
	subaxis(lambda_type);
	ylim(y_lim_vec);
end

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
end






















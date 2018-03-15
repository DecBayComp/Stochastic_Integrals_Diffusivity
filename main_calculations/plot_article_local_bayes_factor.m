%% Plot Bayes factor profile for each simulation and inference model combination


function plot_article_local_bayes_factor(stat_struct, fig_count, bl_save_figures)



%% Constants
load_constants;

% Overrides
marker_size = marker_size + 1;
% line_width = 0;
jitter_amplitude = 0.05;
marker_step = 10;
conventions_count = 4;
line_width = 2;
line_width_marker = 1.25;
legend_font_size = 7;
n_limits = stat_struct.n_limits;
n_limits_count = 3;   % do not plot the full binning
x_lim_vec = [-1.05, 2.05];
output_filename_base = 'K_L';


% y_lim_vec = [-7; 30];

% % % sublabel_x = 0.03;
% % % sublabel_y = 0.08;
% % % x_lim_vec = [x_min, x_max];
% % % y_lim_vec = [-1, 1] * 1.3;
% % % y_lim_vec_profile = [-1, 1] * 0.18;

% % % 
% % % % Overrides
% % % line_width = 1;
% % % 

% Figure size parameters
page_width_frac = 1;
height_factor = 1;

% Subplot parameters
SH = 0.05;
SV = 0.12;
ML = 0.05;
MR = 0.015;
MT = 0.07;
MB = 0.12;
rows = 1;
cols = n_limits_count;


% % % 
% % % % Skip some bins
% % % plot_every = 1;
% % % 
% % % % Label params
% % % sublabel_x = 0.015;
% % % sublabel_y = 1.11;
% % % 
% % % x_tick_increment = 0.2;
% % % 
% % % % Confidence zones
% % % % CONF_LN_K_POSITIVE = [1; 3];
% % % % CONF_LN_K_STRONG = [3; 5];
% % % % CONF_LN_K_VERY_STRONG = [5; 100];
% % % CONF_LN_K_STRONG = 3;
% % % 
% % % % Set confidence zones colors
% % % i=3;
% % % conf_color_sequence = [my_colors(i).White; my_colors(i+1).White; my_colors(i+2).White];
% % % % conf_color_sequence = [233, 243, 216; 234, 246, 255]/255;
% % % 
% % % 
% % % %% Load data
% % % % log_K_L_mean = data_struct.log_K_L_mean;
% % % x_bins_centers = data_struct.x_bins_centers;
% % % bins_count = length(x_bins_centers);
% % % trial_simulation_type = data_struct.trial_simulation_type;
% % % bin_widths = data_struct.x_bins_widths;
% % % bin_borders = [1; 1] * data_struct.x_bins_centers' + [-1/2; 1/2] * bin_widths';
% % % log_K_L_all_trials = data_struct.trials_log_K_L;	% [trials x bins x conventions]
% % % 
% % % % Generate normally distributed random jitter
% % % % jitter = bin_widths * randn(1, conventions_count) * jitter_scale;
% % % jitter = ones(bins_count, 1) * ((0:conventions_count-1) / (conventions_count - 1) - 1/2) * 2 * jitter_scale;
% % % 
% % % 


%% Load data
% mean_log_K_L_avg_neg_x = data_struct.mean_log_K_L_avg_neg_x;
% mean_log_K_L_avg_pos_x = data_struct.mean_log_K_L_avg_pos_x;

% Calculate D' at the first site
% [~, D_prime, ~] = D_func(selected_D_case, x_min, L);
% D_prime_abs = abs(D_prime_ends(1));

% % Calculate alpha for x>0 and x<0 for each convention
% a_ends = f_func(data_struct.selected_f_case, [x_min; x_max], L) / gamma_drag;
% sim_lambda_list = [0, 0.5, 1, 0.5];
% alpha_ends = a_ends * ones(1, length(sim_lambda_list)) + D_prime_ends * sim_lambda_list;
% alpha_ends_over_D = alpha_ends ./ (D_prime_ends * ones(1, 4));

ksi_array = stat_struct.ksi_array;

% Calculate zeta
[D, D_prime] = D_func(selected_D_case, stat_struct.x_bins_centers(stat_struct.middle_bin), L);
b_prime = D_prime / sqrt(2 * D);
zetas = b_prime * sqrt(t_step .* n_limits(1:n_limits_count));
snr = ksi_array' * zetas ./ sqrt(2);

% Create jitter to distinguish curves
jitter_array = (-1.5:1.5) * jitter_amplitude;

% Determine y limits
y_lim_vec(1) = min(reshape(stat_struct.log_K_L_mean(:, 1:n_limits_count, 1:conventions_count) - stat_struct.log_K_L_eb(:, 1:n_limits_count, 1:conventions_count), [], 1));
y_lim_vec(2) = max(reshape(stat_struct.log_K_L_mean(:, 1:n_limits_count, 1:conventions_count) + stat_struct.log_K_L_eb(:, 1:n_limits_count, 1:conventions_count), [], 1));



%% Plot
% Initalize plot
h_fig = figure(fig_count);
set_article_figure_size(h_fig, rows, page_width_frac, height_factor);
clf;


% Initalize subplots
subaxis(rows, cols, 1, 'SH', SH, 'SV', SV, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);

for lim_ind = 1:n_limits_count
    
    subaxis(n_limits_count - lim_ind + 1);
    hold on;
    str_legend = cell(conventions_count, 1);
    hand_curves = zeros(conventions_count,1);
    for convention = 1:conventions_count
%         plot(ksi_array + jitter_array(convention), stat_struct.log_K_L_mean(:, lim_ind, convention), ...
%             strcat('-'), 'color', color_sequence(convention, :),  'LineWidth', line_width, 'markers', marker_size);

        % Put markers
        hand_curves(convention) = plot(ksi_array(1:marker_step:end), stat_struct.log_K_L_mean(1:marker_step:end, lim_ind, convention), ...
             strcat(markers_list{convention}), 'color', color_sequence(convention, :),  'LineWidth', line_width_marker, 'markers', marker_size);

        % Load data
        x = ksi_array;
        avg = stat_struct.log_K_L_mean(:, lim_ind, convention);
        err = stat_struct.log_K_L_eb(:, lim_ind, convention)';
        curve = shadedErrorBar(x, avg, err, 'lineprops', {'color', color_sequence(convention, :), 'LineWidth', line_width}, 'patchSaturation', 0.125);
        %'color', color_sequence(convention, :),  'LineWidth', line_width_error, 'markers', marker_size, markers_list{convention});
        
        % Remove edge
        set(curve.edge, 'visible', 'off');
        
        str_legend{convention} = conventions_names{convention};
%         hand_curves(convention) = curve.mainLine;
    end

    % % Adjust
    % x_lim_vec = xlim();
    % x_lim_vec(1) = 0;
    xlim(x_lim_vec);
    ylim(y_lim_vec);
    set(gca, 'FontSize', font_size);

    xlabel('$\alpha / bb''$', 'Interpreter', 'latex', 'FontSize', font_size);
    ylabel('$\langle \ln K_L \rangle$', 'Interpreter', 'latex', 'FontSize', font_size);
    title(sprintf('$n = %i$', n_limits(lim_ind)), 'Interpreter', 'latex', 'FontSize', subplot_label_font_size);
    
    % Theory
    h_theor = plot(xlim(), [0, 0], 'LineWidth', line_width_theor, 'color', axes_color);
    uistack(h_theor, 'bottom');
    
    % Legend
	if lim_ind == 3
		legend(hand_curves, str_legend, 'location', 'north', 'FontSize', legend_font_size);
    end
end






















% % % 
% % % %% Plot
% % % % Initalize plot
% % % h_fig = figure(fig_count);
% % % set_article_figure_size(h_fig, rows, page_width_frac, height_factor);
% % % clf;
% % % 
% % % % Initalize subplots
% % % subaxis(rows, cols, 1, 'SH', SH, 'SV', SV, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);
% % % % for lambda_type = 1:lambda_types_count
% % % 
% % % 
% % % 
% % % %% Calculate
% % % mean_log_K_L = zeros(lambda_types_count, bins_count, conventions_count);
% % % std_log_K_L = zeros(lambda_types_count, bins_count, conventions_count);
% % % eb_log_K_L = zeros(lambda_types_count, bins_count, conventions_count);	% error bars
% % % mean_evidence = zeros(lambda_types_count, bins_count, conventions_count);
% % % for lambda_type = 1:lambda_types_count
% % % 	% Extract current-lambda-type results
% % % 	log_K_L = log_K_L_all_trials(trial_simulation_type == lambda_type, :, :);
% % % 	
% % % 	% Calculate mean and std of the log(K)! over all trials of the same type
% % % 	mean_log_K_L(lambda_type, :, :) = mean(log_K_L, 1);
% % % 	std_log_K_L(lambda_type, :, :) = std(log_K_L, [], 1);
% % % 	eb_log_K_L(lambda_type, :, :) = std_log_K_L(lambda_type, :, :) * sqrt(2) * erfinv(0.95);
% % % end
% % % 
% % % % Evaluate evidence for the presence of force on 3-values scale
% % % mean_evidence(mean_log_K_L >= CONF_LN_K_STRONG) = 1;
% % % mean_evidence(mean_log_K_L <= -CONF_LN_K_STRONG) = -1;
% % % 
% % % 
% % % 	
% % % for lambda_type = 1:lambda_types_count	
% % % 	%% Plot
% % % 	subaxis(lambda_type);
% % % 	hold on;
% % %     
% % %     % Theory
% % % 	h_theor = plot(x_lim_vec, [0,0], '-', 'linewidth', line_width_theor, 'color', axes_color);
% % %     
% % % 	str_legend = {};
% % %     conv_plots = zeros(1, conventions_count);
% % %     for convention = 1:conventions_count
% % % 		conv_plots(convention) = plot(x_bins_centers, squeeze(mean_evidence(lambda_type, :, convention)) + jitter(:, convention)', strcat('', markers_list{convention}),...
% % % 			'color', color_sequence(convention, :), 'markers', marker_size, 'linewidth', line_width);
% % % 		str_legend{length(str_legend) + 1} = conventions_names{convention};
% % %     end
% % % % 	x_lim_vec = xlim();
% % % 	xlim(x_lim_vec);
% % % 
% % % 	
% % % 
% % % 	% Adjust
% % % 	ylim(y_lim_vec);
% % % % 	ylim('auto');
% % % 
% % % 	% Subplot label
% % % 	text(sublabel_x, sublabel_y, char('A' + lambda_type - 1 + 4 * bl_force), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', subplot_label_font_size);
% % % 
% % % 	% Modify ticks
% % % 	set(gca,'xtick', x_min:x_tick_increment:x_max);
% % % 	set(gca,'ytick', [])
% % % 	if lambda_type == 1
% % % 		set(gca,'ytick', [-1, 0, 1]);
% % %     end
% % % 
% % % 
% % % 	% Axes labels
% % % 	xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
% % % 	if lambda_type == 1 && ~bl_force
% % % 		ylabel({'No force', '$\langle \ln K_L \rangle$'}, 'interpreter', 'latex');
% % % 	elseif lambda_type == 1 && bl_force
% % % 		ylabel({'Force', '$\langle \ln K_L \rangle$'}, 'interpreter', 'latex');
% % % 	end
% % % 
% % % 	% Title
% % % 	str_title = sprintf('%s sim.', lambda_types_tex_names{lambda_type});
% % % 	title(str_title, 'interpreter', 'latex');
% % % 	
% % % 	% Legend
% % % 	if lambda_type == 2 && ~data_struct.bl_force
% % % 		legend(conv_plots, str_legend, 'location', 'northeast', 'fontsize', legend_font_size);
% % % 	
% % % 	elseif lambda_type == 2 && data_struct.bl_force
% % % 		legend(conv_plots, str_legend, 'location', 'southwest', 'fontsize', legend_font_size);
% % %     end
% % % 	
% % % 	% Add confidence zones
% % % 	x_lim_vec = xlim();
% % % 	
% % % 	% Color bin borders
% % % 	color_bins(bin_borders, y_lim_vec, bin_color);
% % % 	
% % % 	% Put axes on top
% % % 	set(gca, 'Layer', 'top');
% % % 
% % % end
% % % 
% % % 
% % % 


%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize', [fig_pos(3), fig_pos(4)]);
% Set filename
output_filename = strcat(output_filename_base, '.pdf');
output_full_path = strcat(output_figures_folder, output_filename);
if bl_save_figures
    print(h_fig, output_full_path, '-dpdf', '-r0');
end










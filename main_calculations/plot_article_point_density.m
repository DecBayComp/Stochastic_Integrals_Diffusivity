%% Plot the distribution of point numbers in bins


function plot_article_point_density(stat_struct, fig_count, bl_save_figures)


%% Constants
load_constants;
x_tick_increment = 0.5;
y_lim_vec = [0, 0.25];
output_filename_base = 'point_density';

% Figure size parameters
page_width_frac = 0.25;
height_factor = 0.6;

% Label params
sublabel_x = 0.0;
sublabel_y = 1.23;

% Define plot colors
load_color_scheme;
color_sequence = [standard_colors(1).DeepBlue; my_colors(5).Green; my_colors(1).Orange; my_colors(1).WarmBrown];
% bin_color = [standard_colors(1).LightBlue, transp];



%% Initialize figure
h_fig = figure(fig_count);
set_article_figure_size(h_fig, 1, page_width_frac, height_factor);
% tightfig(h_fig);
clf;
hold on;



%% Initialize
bins_number = stat_struct.x_bins_number;
bin_widths = stat_struct.x_bins_widths;
bin_borders = [1; 1] * stat_struct.x_bins_centers' + [-1/2; 1/2] * bin_widths';

ksi_count = length(stat_struct.ksi_array);
n_limits_count = length(n_limits);

% b_profile = stat_struct.b_theor_fine_data;

x_lim_vec = [x_min, x_max];



%% Plot point number in bins n_j
norm_points_density = zeros(ksi_count, bins_number);
lim_ind = n_limits_count;
str_legend = {};
for ksi_ind = 1:ksi_count
	% Load data
	norm_points_density(ksi_ind, :) = stat_struct.n_j_mean(ksi_ind, lim_ind, :);
	
	% Normalize to bin width and to 1
	norm_points_density(ksi_ind, :) = norm_points_density(ksi_ind, :) ./ bin_widths';
	norm_points_density(ksi_ind, :) = norm_points_density(ksi_ind, :) ./ sum(norm_points_density(ksi_ind, :));
	
	% Plot
    plot(stat_struct.x_bins_centers,  norm_points_density(ksi_ind, :), 'LineWidth', line_width - 0.8);
%         strcat('-', markers_list{ksi_ind}), 'color', color_sequence(ksi_ind, :),  , 'markers', marker_size);
%     str_legend{end + 1} = lambda_types_names{ksi_ind};
end

% Adjust plot
xlim(x_lim_vec);
box on;
xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
ylabel('Point density', 'interpreter', 'latex');
%  
% if ~bl_force
% 	title('No force', 'interpreter', 'latex');
% else
% % 	set(gca,'YTickLabel',[]);
% 	title('Force', 'interpreter', 'latex');
% end

% Modify ticks
set(gca, 'FontSize', font_size);
set(gca,'xtick', x_min:x_tick_increment:x_max);

% Legend 
% legend(str_legend, 'Location', 'southwest', 'FontSize', legend_font_size);

% Subplot label
% if ~bl_force
% 	chr_label = 'A';
% else
% 	chr_label = 'B';
% end
% text(sublabel_x, sublabel_y, chr_label, 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', subplot_label_font_size);

% Adjust y limits
% y_lim_vec = ylim();
% y_lim_vec(1) = 0;
ylim(y_lim_vec);



% % Normalize and plot simulated diffusivity profile
% b_profile = b_profile - min(b_profile);
% b_profile = b_profile / max(b_profile) * y_lim_vec(2);
% h_theor_center = plot(data_struct.x_fine_mesh, b_profile, '-k', 'LineWidth', line_width_theor);
% % Send curve back
% uistack(h_theor_center, 'bottom');

% Color bin borders
color_bins(bin_borders, y_lim_vec, bin_color);



% %% Save figure
% % Prepare printer
% h_fig.PaperPositionMode = 'auto';
% h_fig.Units = 'Inches';
% fig_pos = h_fig.Position;
% set(h_fig, 'PaperUnits','Inches','PaperSize', [fig_pos(3), fig_pos(4)]);
% 
% % Set filename
% output_filename = strcat(output_filename_base, '_', data_struct.str_force, '.pdf');
% output_full_path = strcat(output_figures_folder, output_filename);
% if bl_save_figures
%     print(h_fig, output_full_path, '-dpdf', '-r0');
% end




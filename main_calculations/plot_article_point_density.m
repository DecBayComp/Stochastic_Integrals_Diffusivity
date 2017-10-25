%% Plot the distribution of point numbers in bins


function plot_article_point_density(data_struct, fig_count, bl_save_figures)


%% Constants
load_constants;
x_tick_increment = 0.1;
output_filename = 'Point_density.pdf';

% Define plot colors
load_color_scheme;
color_sequence = [standard_colors(1).DeepBlue; my_colors(5).Green; my_colors(1).Orange; my_colors(1).WarmBrown];
% bin_color = [standard_colors(1).LightBlue, transp];



%% Initialize figure
h_fig = figure(fig_count);
set_article_figure_size(h_fig, 1, 1, 1);
clf;
hold on;



%% Initialize
bins_number = data_struct.x_bins_number;
bin_widths = data_struct.x_bins_widths;
bin_borders = [1; 1] * data_struct.x_bins_centers' + [-1/2; 1/2] * bin_widths';

b_profile = data_struct.b_theor_fine_data;

x_lim_vec = [x_min, x_max];



%% Plot point number in bins n_j
norm_points_density = zeros(lambda_types_count, bins_number);
str_legend = {};
for lambda_type = 1:lambda_types_count
	% Load data
	norm_points_density(lambda_type, :) = data_struct.n_j_mean(lambda_type, :);
	
	% Normalize to total points number
	norm_points_density(lambda_type, :) = norm_points_density(lambda_type, :) / sum(norm_points_density(lambda_type, :), 2) * 100;
	
	% Plot
    plot(data_struct.x_bins_centers,  norm_points_density(lambda_type, :),...
        strcat('-', markers_list{lambda_type}), 'color', color_sequence(lambda_type, :),  'LineWidth', line_width, 'markers', marker_size);
    str_legend{end + 1} = lambda_types_names{lambda_type};
end;

% Adjust plot
xlim(x_lim_vec);
box on;
xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
ylabel('Points in bin, \%', 'interpreter', 'latex');
title('Point density and bin locations', 'interpreter', 'latex');

% Modify ticks
set(gca,'xtick', x_min:x_tick_increment:x_max);

% Legend 
legend(str_legend, 'Location', 'northwest');

% Normalize and plot simulated diffusivity profile
y_lim_vec = ylim();
b_profile = b_profile - min(b_profile);
b_profile = b_profile / max(b_profile) * y_lim_vec(2);
h_theor_center = plot(data_struct.x_fine_mesh, b_profile, '-k', 'LineWidth', line_width_theor);
% Send curve back
uistack(h_theor_center, 'bottom');

% Color bin borders
color_bins(bin_borders, y_lim_vec, bin_color);



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







function plot_article_mean_a_profile_and_fail_rate(data_struct, trials_data, fig_count, bl_save_figures)
%% === Plot mean a profile for each inference type ===


%% Constants
load_constants;
sublabel_x = 0.03;
sublabel_y = 0.08;
x_lim_vec = [x_min, x_max];
y_lim_vec_FR = [-2, 102];
y_lim_vec_profile = [-1, 1] * 0.07;
output_filename_base = 'a';
% Subplot parameters
SH = 0.06;
SV = 0.1;
ML = 0.07;
MR = 0.01;
MT = 0.1;
MB = 0.15;
rows = 1;
cols = 2;

% Skip some bins
plot_every = 1;

% Label params
sublabel_x = 0.015;
sublabel_y = 1.1;

x_tick_increment = 0.2;

selected_x_over_L = -0.25;



%% Plot
% Initalize plot
h_fig = figure(fig_count);
set_article_figure_size(h_fig, rows, 2, 1);
clf;
% Initalize subplots
subaxis(rows, cols, 1, 'SH', SH, 'SV', SV, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);
% for lambda_type = 1:lambda_types_count
lambda_type = enum_lambda_Hanggi;
    


%% == (1): Profile plot ==
subaxis(1);
hold on;
% Plot each convention
for convention = 1:conventions_count
	plot(data_struct.x_bins_centers(1:plot_every:end),  data_struct.MAP_a_mean(lambda_type, 1:plot_every:end, convention, 1),...
		strcat('-', markers_list{convention}), 'color', color_sequence(convention, :), 'LineWidth', line_width, 'markers', marker_size);
end;
% True profile (theory)
h_theor = plot(data_struct.x_fine_mesh, data_struct.a_theor_fine_data, '--k', 'LineWidth', line_width);

%% Adjust
xlim(x_lim_vec);
ylim(y_lim_vec_profile);
box on;

xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
%     if lambda_type == 1
	ylabel('$\langle \hat a \rangle, \mu \mathrm{m/s}$', 'interpreter', 'latex');
%     end;
title('Mean local force profile', 'interpreter', 'latex');
%     str_title = {'$\lambda^* = 0$', '$\lambda^* = 0.5$', '$\lambda^* = 1$', 'Random $\lambda^*$'};
%     title(str_title{lambda_type}, 'interpreter', 'latex');

% Subplot label
text(sublabel_x, sublabel_y, 'A', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', subplot_label_font_size);

% Push theoretical curve back
uistack(h_theor, 'bottom');

% Modify ticks
set(gca,'xtick', x_min:x_tick_increment:x_max);



%% == (2): a profile in bin ==
lambda_type = enum_lambda_Hanggi;

subaxis(2);
plot_article_a_profile_in_bin(data_struct, trials_data, lambda_type, selected_x_over_L);

% Subplot label
text(sublabel_x, sublabel_y, 'B', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', subplot_label_font_size);
	


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










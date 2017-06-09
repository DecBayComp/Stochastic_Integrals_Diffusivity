

function plot_article_D(data_struct, trials_data, fig_count, bl_save_figures)


%% Constants
load_constants;
bin_plot_step = 3;
rows = 1;
cols = 3;
x_lim_vec = [0, x_max];
% Subplots params
SH = 0.07;
SV = 0.125;
ML = 0.07;
MR = 0.015;
MT = 0.09;
MB = 0.18;
% Label params
sublabel_x = 0.9;
sublabel_y = 0.135;
output_filename = 'D.pdf';


%% == (A): Plot diffusivity profile ==
% Initialize
h_fig = figure(fig_count);
set_article_figure_size(h_fig, rows, 2, 1);
clf;
h_sub = subaxis(rows, cols, 1, 'SH', SH, 'SV', SV, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);
hold on;
% Plot
str_legend = {};
for lambda_type = 1:lambda_types_count
    plot(data_struct.x_bins_centers(1:bin_plot_step:end),  data_struct.MAP_D_mean(lambda_type, 1:bin_plot_step:end, 1),...
        strcat('-', markers_list{lambda_type}), 'color', color_sequence(lambda_type, :),  'LineWidth', line_width, 'markers', marker_size);
    str_legend{end + 1} = lambda_types_names{lambda_type};
end;
% Theory
h_theor = plot(data_struct.x_fine_mesh, data_struct.D_theor_fine_data, '--k', 'LineWidth', line_width);
% Adjust
xlim(x_lim_vec);
y_lim_vec = [0.004, 0.016];
ylim(y_lim_vec);
box on;
xlabel('$x$', 'interpreter', 'latex');
ylabel('$D$', 'interpreter', 'latex');
title('Average $D$ profile', 'interpreter', 'latex');
% Subplot label
text(sublabel_x, sublabel_y, '(a)', 'Units', 'Normalized', 'VerticalAlignment', 'Top');
% Legend
% str_legend = lambda_types_names;
h_leg = legend(str_legend, 'location', 'northwest', 'FontSize', legend_font_size);
legend boxon;
% Send theoretical curve back
uistack(h_theor, 'bottom');


%% == (B): Fail rate ==
subaxis(rows, cols, 2);
hold on;
y_lim_vec = [0, 1];
str_legend = {};
for lambda_type = 1:lambda_types_count
    plot(data_struct.x_bins_centers,  data_struct.UR_D(lambda_type, :),...
        strcat('-', markers_list{lambda_type}), 'color', color_sequence(lambda_type, :),  'LineWidth', line_width, 'markers', marker_size);
    str_legend{end + 1} = lambda_types_names{lambda_type};
end;
% Plot the used confidence level
h_conf = plot(x_lim_vec, [1, 1] * (1 - CONF_LEVEL), 'k--', 'linewidth', line_width);
% Adjust
xlim(x_lim_vec);
ylim(y_lim_vec);
box on;
xlabel('$x$', 'interpreter', 'latex');
ylabel('Fail rate', 'interpreter', 'latex');
title('Average fail rate', 'interpreter', 'latex');
% Legend
% str_legend = lambda_types_names;
legend(str_legend, 'location', 'northwest', 'FontSize', legend_font_size);
legend boxon;
% Subplot label
text(sublabel_x, sublabel_y, '(b)', 'Units', 'Normalized', 'VerticalAlignment', 'Top');
% Send confidence level back
uistack(h_conf, 'bottom');


%% == (C): D' profile ==
%% Initialize
subaxis(rows, cols, 3);
hold on;
y_lim_vec = [-1, 0.75] * 1;
%% Calculate
% Calculate simple D' with MAP D using a finite elements scheme
% Choose one lambda and corresponding data_struct
tmp_data_struct = trials_data{data_struct.trial_first_simulation_type_index(enum_lambda_Stratonovich)};
x_bins_steps = tmp_data_struct.x_bins_centers(2:end) - tmp_data_struct.x_bins_centers(1:end - 1);
simple_D_grad = (tmp_data_struct.MAP_D(2:end, 1) - tmp_data_struct.MAP_D(1:end-1, 1)) ./ x_bins_steps;
x_grad_mesh = tmp_data_struct.x_grad_mesh;
%% Plot
% Simple difference
plot(x_grad_mesh, simple_D_grad, markers_list{1}, 'color', color_sequence(1, :), 'LineWidth', line_width, 'markers', marker_size);
% Regularized gradient
plot(x_grad_mesh, tmp_data_struct.MAP_D_grad_regular,  markers_list{2}, 'color', color_sequence(2, :), 'LineWidth', line_width, 'markers', marker_size);
% Regularized interpolated gradient
plot(tmp_data_struct.x_bins_centers, tmp_data_struct.MAP_D_grad_regular_interp, 'color', color_sequence(3, :), 'LineWidth', line_width, 'markers', marker_size);
% Theory
h_theor = plot(tmp_data_struct.x_fine_mesh, tmp_data_struct.D_grad_theor_fine_data, 'k--', 'LineWidth', line_width);
% Adjust
xlabel('$x$', 'interpreter', 'latex');
ylabel('$\nabla D$', 'interpreter', 'latex');
xlim(x_lim_vec);
ylim(y_lim_vec);
box on;
title(sprintf('$\\nabla D$ profile for $\\lambda^* = %.2f$', tmp_data_struct.lambda), 'interpreter', 'latex');
% Sublabel
text(sublabel_x, sublabel_y, '(c)', 'Units', 'Normalized', 'VerticalAlignment', 'Top');
% Legend
str_legend_local = {'FD', 'R', 'RI'};
legend(str_legend_local, 'location', 'southwest', 'interpreter', 'latex', 'FontSize', legend_font_size);
% Send theoretical curve back
uistack(h_theor, 'bottom');


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












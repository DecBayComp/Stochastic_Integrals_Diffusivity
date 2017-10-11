

function plot_article_b(data_struct, trials_data, fig_count, bl_save_figures)


%% Constants
load_constants;
bin_plot_step = 1;	% 3
rows = 2;
cols = 2;
x_lim_vec = [x_min, x_max];
w = 10.0;
% Subplots params
SH = 0.07;
SV = 0.125;
ML = 0.08;
MR = 0.0125;
MT = 0.04;
MB = 0.075;
% Label params
sublabel_x = 0.025;
sublabel_y = 0.07;
output_filename = 'b.pdf';


%% Calculate
% Filter indices from one best period
x_left = (1/2 + 2*1)/w;
x_right = (1/2 + 2*2)/w;
indices = data_struct.x_bins_centers >= x_left & data_struct.x_bins_centers <= x_right;


%% == (A): Plot diffusivity profile ==
%% Calculate the mean theoretical values in bins
bins_borders = [1; 1] * data_struct.x_bins_centers' + [-1/2; 1/2] * data_struct.x_bins_widths';
% Calculate the theoretical anti-derivative of D at the borders
[~, ~, ~, D_antider_theor] = D_func(selected_D_case, bins_borders, L);
% Calculate theor mean D and b in bin
D_mean_theor = (D_antider_theor(2, :) - D_antider_theor(1, :)) ./ (bins_borders(2, :) - bins_borders(1, :));
b_mean_theor = sqrt(2 * D_mean_theor);

% Initialize
y_lim_vec_A = [0.09, 0.18];
h_fig = figure(fig_count);
set_article_figure_size(h_fig, rows, 2, 1);
clf;
h_sub = subaxis(rows, cols, 1, 'SH', SH, 'SV', SV, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);
hold on;
% Plot
str_legend = {};
for lambda_type = 1:lambda_types_count
    plot(data_struct.x_bins_centers(1:bin_plot_step:end),  data_struct.MAP_b_mean(lambda_type, 1:bin_plot_step:end, 1),...
        strcat('-', markers_list{lambda_type}), 'color', color_sequence(lambda_type, :),  'LineWidth', line_width, 'markers', marker_size);
    str_legend{end + 1} = lambda_types_names{lambda_type};
end;
%% Theory
% Values in center
h_theor_center = plot(data_struct.x_fine_mesh, data_struct.b_theor_fine_data, '-k', 'LineWidth', line_width);
% Mean values in bins
plot(data_struct.x_bins_centers, b_mean_theor, '--k', 'LineWidth', line_width);

% Adjust
xlim(x_lim_vec);
ylim(y_lim_vec_A);
box on;
xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
ylabel('$\langle b \rangle$, $\mu\mathrm{m \cdot s^{-1/2}}$', 'interpreter', 'latex');
title('Average diffusivity profile', 'interpreter', 'latex');
% Subplot label
text(sublabel_x, sublabel_y, '(a)', 'Units', 'Normalized', 'VerticalAlignment', 'Top');
% Legend
% str_legend = lambda_types_names;
h_leg = legend(str_legend, 'location', 'south', 'FontSize', legend_font_size);
legend boxon;
grid on;
% Send theoretical curve back
uistack(h_theor_center, 'bottom');


%% == (B): Fail rate ==
subaxis(rows, cols, 2);
hold on;
y_lim_vec = [-2, 102];
str_legend = {};
for lambda_type = 1:lambda_types_count
    plot(data_struct.x_bins_centers,  data_struct.UR_b(lambda_type, :) * 100,...
        strcat('-', markers_list{lambda_type}), 'color', color_sequence(lambda_type, :),  'LineWidth', line_width, 'markers', marker_size);
    str_legend{end + 1} = lambda_types_names{lambda_type};
end;
% Plot the used confidence level
h_conf = plot(x_lim_vec, [1, 1] * (1 - CONF_LEVEL) * 100, 'k--', 'linewidth', line_width);
% Adjust
xlim(x_lim_vec);
ylim(y_lim_vec);
box on;
xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
ylabel('Fail rate, \%', 'interpreter', 'latex');
title('Fail rate', 'interpreter', 'latex');
% Legend
% str_legend = lambda_types_names;
legend(str_legend, 'location', 'north', 'FontSize', legend_font_size);
legend boxon;
% Subplot label
text(sublabel_x, sublabel_y, '(b)', 'Units', 'Normalized', 'VerticalAlignment', 'Top');
% Send confidence level back
uistack(h_conf, 'bottom');



%% == (C): b bias ==

%% Calculate the theoretically expected bias based on integral series
% Load values
x_mesh = data_struct.x_bins_centers;
bin_sizes = data_struct.x_bins_widths;
b_theor_data = data_struct.b_theor_data;
mean_jumps = data_struct.mean_jump_bins_all_trials';

% b' part
lambdas = [0, 0.5, 1];
b_mean_series = 2 * b_theor_data(:, 2).^2 .* (bin_sizes / 2).^2 ./ 6 ./ b_theor_data(:, 1);
b_mean_series = b_mean_series * (2 * lambdas - 2);

% b'' part
b_mean_series = b_mean_series + b_theor_data(:, 3) .* (bin_sizes / 2).^2 / 6 * [1, 1, 1];

% Add the central value and subtract the theoretical mean
b_theor_bias = b_mean_series + (data_struct.b_theor_data(:, 1) - b_mean_theor') * [1, 1, 1];

%% Initialize
subaxis(rows, cols, 3);
hold on;
% x_lim_vec_C = [-0.19,0.16];
y_lim_vec = [-1, 1] * 7.5e-3;
str_legend = {};

%% Plot
for lambda_type = 1:lambda_types_count
    plot(data_struct.x_bins_centers,...
        (data_struct.MAP_b_mean(lambda_type, :, 1) - b_mean_theor), strcat('-', markers_list{lambda_type}),...
        'markers', marker_size, 'LineWidth', line_width, 'color', color_sequence(lambda_type, :));
    str_legend{end + 1} = lambda_types_names{lambda_type};
end;

%% Legend
legend(str_legend, 'location', 'northwest', 'FontSize', legend_font_size);
legend boxon;

%% Theory
for lambda_ind = 1:length(lambdas)
	plot(data_struct.x_bins_centers, b_theor_bias(:, lambda_ind), 'k--', 'color', color_sequence(lambda_ind, :));
end;
%h_conf = plot(x_lim_vec, [1, 1] * (1 - CONF_LEVEL) * 100, 'k--', 'linewidth', line_width);

%% Plot zero bias
h_theor_0 = plot(x_lim_vec, 0 * x_lim_vec, 'k--');

%% Adjust
xlim(x_lim_vec);
ylim(y_lim_vec);
xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
ylabel('$\langle \delta b \rangle$, $\mu\mathrm{m \cdot s^{-1/2}}$', 'interpreter', 'latex');
text(sublabel_x, sublabel_y, strcat('(c)'), 'Units', 'Normalized', 'VerticalAlignment', 'Top');
title('Average diffusivity bias', 'interpreter', 'latex');
grid on;
% Reorder curves
uistack([h_theor_0], 'bottom');



%% == (D): D' profile ==
%% Initialize
subaxis(rows, cols, 4);
hold on;
y_lim_vec = [-1, 1] * 0.23;
%% Calculate
% Calculate simple bb' with MAP b using a finite elements scheme
% Choose one lambda and corresponding data_struct
tmp_data_struct = trials_data{data_struct.trial_first_simulation_type_index(enum_lambda_Stratonovich)};
x_bins_steps = tmp_data_struct.x_bins_centers(2:end) - tmp_data_struct.x_bins_centers(1:end - 1);
b_squared_over_2 = tmp_data_struct.MAP_b.^2 / 2;
simple_bb_prime = (b_squared_over_2(2:end, 1) - b_squared_over_2(1:end-1, 1)) ./ x_bins_steps;
x_grad_mesh = tmp_data_struct.x_grad_mesh;
%% Plot
% Simple difference
plot(x_grad_mesh, simple_bb_prime, markers_list{1}, 'color', color_sequence(1, :), 'LineWidth', line_width, 'markers', marker_size);
% Regularized gradient
plot(x_grad_mesh, tmp_data_struct.MAP_bb_prime_regular,  markers_list{2}, 'color', color_sequence(2, :), 'LineWidth', line_width, 'markers', marker_size);
% Regularized interpolated gradient
plot(tmp_data_struct.x_bins_centers, tmp_data_struct.MAP_bb_prime_regular_interp, 'color', color_sequence(3, :), 'LineWidth', line_width, 'markers', marker_size);
% Theory
h_theor = plot(tmp_data_struct.x_fine_mesh, tmp_data_struct.bb_prime_theor_fine_data, 'k--', 'LineWidth', line_width);
% Adjust
xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
ylabel('$bb''$, $\mu \mathrm{m/s}$', 'interpreter', 'latex');
xlim(x_lim_vec);
ylim(y_lim_vec);
box on;
grid on;
title(sprintf('$bb''$ profile for $\\lambda^* = %.2f$', tmp_data_struct.lambda), 'interpreter', 'latex');
% Sublabel
text(sublabel_x, sublabel_y, '(d)', 'Units', 'Normalized', 'VerticalAlignment', 'Top');
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












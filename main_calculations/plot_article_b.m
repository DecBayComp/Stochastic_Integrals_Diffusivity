

function plot_article_b(data_struct, trials_data, fig_count, bl_save_figures)


%% Constants
load_constants;
bin_plot_step = 3;
rows = 2;
cols = 2;
x_lim_vec = [0, x_max];
w = 10.0;
% Subplots params
SH = 0.07;
SV = 0.125;
ML = 0.07;
MR = 0.0125;
MT = 0.04;
MB = 0.075;
% Label params
sublabel_x = 0.025;
sublabel_y = 0.07;
output_filename = 'D.pdf';


%% Calculate
% Filter indices from one best period
x_left = (1/2 + 2*1)/w;
x_right = (1/2 + 2*2)/w;
indices = data_struct.x_bins_centers >= x_left & data_struct.x_bins_centers <= x_right;


%% == (A): Plot diffusivity profile ==
% Initialize
y_lim_vec = [0.004, 0.01525];
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
ylim(y_lim_vec);
box on;
xlabel('$x$', 'interpreter', 'latex');
ylabel('$<D>$', 'interpreter', 'latex');
title('Average $D$ profile', 'interpreter', 'latex');
% Subplot label
text(sublabel_x, sublabel_y, '(a)', 'Units', 'Normalized', 'VerticalAlignment', 'Top');
% Legend
% str_legend = lambda_types_names;
h_leg = legend(str_legend, 'location', 'south', 'FontSize', legend_font_size);
legend boxon;
% Send theoretical curve back
uistack(h_theor, 'bottom');


%% == (B): Fail rate ==
subaxis(rows, cols, 2);
hold on;
y_lim_vec = [-2, 102];
str_legend = {};
for lambda_type = 1:lambda_types_count
    plot(data_struct.x_bins_centers,  data_struct.UR_D(lambda_type, :) * 100,...
        strcat('-', markers_list{lambda_type}), 'color', color_sequence(lambda_type, :),  'LineWidth', line_width, 'markers', marker_size);
    str_legend{end + 1} = lambda_types_names{lambda_type};
end;
% Plot the used confidence level
h_conf = plot(x_lim_vec, [1, 1] * (1 - CONF_LEVEL) * 100, 'k--', 'linewidth', line_width);
% Adjust
xlim(x_lim_vec);
ylim(y_lim_vec);
box on;
xlabel('$x$', 'interpreter', 'latex');
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


%% == (C): D bias ==
%% Initialize
subaxis(rows, cols, 3);
hold on;
x_lim_vec_C = [-0.19,0.16];
y_lim_vec = [-1.1, 0.65] * 1e-3;
str_legend = {};
%% Plot
for lambda_type = 1:lambda_types_count
    plot(data_struct.MAP_D_grad_regular_interp_mean(lambda_type, indices) * kBT,...
        (data_struct.MAP_D_mean(lambda_type, indices, 1) - data_struct.D_theor_data(indices)'), markers_list{lambda_type},...
        'markers', marker_size, 'LineWidth', line_width, 'color', color_sequence(lambda_type, :));
    str_legend{end + 1} = lambda_types_names{lambda_type};
end;
%% Legend
legend(str_legend, 'location', 'southwest', 'FontSize', legend_font_size);
legend boxon;
%% Theory
xlim(x_lim_vec_C);
ylim(y_lim_vec);
% y = 0
h_theor_0 = plot(x_lim_vec_C, 0 * x_lim_vec_C, 'k--');
%% Adjust
xlabel('$k_\mathrm{B}T<\nabla D>$', 'interpreter', 'latex');
ylabel('$<\delta D>$', 'interpreter', 'latex');
text(sublabel_x, sublabel_y, strcat('(c)'), 'Units', 'Normalized', 'VerticalAlignment', 'Top');
title('Average $D$ bias', 'interpreter', 'latex');
% Reorder curves
uistack([h_theor_0], 'bottom');


%% == (D): D' profile ==
%% Initialize
subaxis(rows, cols, 4);
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
text(sublabel_x, sublabel_y, '(d)', 'Units', 'Normalized', 'VerticalAlignment', 'Top');
% Legend
str_legend_local = {'FD', 'R', 'RI'};
legend(str_legend_local, 'location', 'northwest', 'interpreter', 'latex', 'FontSize', legend_font_size);
% Send theoretical curve back
uistack(h_theor, 'bottom');



% Plot bias in D against the local d gradient
% figure(9);


%     x_lim_vec = xlim();
% 
%     %% Identity-like lines (theory)
%     % y = x
%     h_theor_1 = plot(x_lim_vec, x_lim_vec, 'k--');
%     % y = x/2
%     h_theor_2 = plot(x_lim_vec, 1/2 * x_lim_vec, 'k--');
%     % y = 0
%     h_theor_0 = plot(x_lim_vec, 0 * x_lim_vec, 'k--');
%     % y = -x/2
%     h_theor_m2 = plot(x_lim_vec, -1/2 * x_lim_vec, 'k--');
%     % y = -x
%     h_theor_m2 = plot(x_lim_vec, -1 * x_lim_vec, 'k--');
    

%     % Legend
%     if lambda_type == 1
%         legend(str_legend, 'location', 'northeast', 'FontSize', legend_font_size);
%     end;
    
%     xlim(x_lim_vec);
%     ylim(y_lim_vec);
    
    % Labels
   
%     title_str = {'$\lambda^* = 0$', '$\lambda^* = 0.5$', '$\lambda^* = 1$', 'Random $\lambda^*$'};
%     title(title_str{lambda_type}, 'interpreter', 'latex');
    % Subplot label
    
    
    % Reorder curves
%     uistack([h_theor_0, h_theor_m2, h_theor_2, h_theor_1], 'bottom');




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












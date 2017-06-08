


function plot_article_fD_bias(data_struct, fig_count, bl_save_figures)
%% === Plot fD bias as a function fo the gradient ===


%% Constants
load_constants;
sublabel_x = 0.02;
sublabel_y = 0.99;
x_lim_vec = [-1, 1] * 0.17;
y_lim_vec = [-1, 1] * 0.18;
output_filename = 'Force_bias_vs_gradient.pdf';
% Subplots parameters
spacing = 0.08;
ML = 0.07;
MR = 0.02;
MT = 0.08;
MB = 0.15;


%% Calculate
% Filter indices from one best period
w = 10.0;
x_left = (1/2 + 2*1)/w;
x_right = (1/2 + 2*2)/w;
indices = data_struct.x_bins_centers >= x_left & data_struct.x_bins_centers <= x_right;

    
%% Plot
h_fig = figure(fig_count);
set_article_figure_size(h_fig, 1, 2, 1);
clf;
subaxis(1, lambda_types_count, 1, 'Spacing', spacing, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);
for lambda_type = 1:lambda_types_count
    subaxis(1, lambda_types_count, lambda_type);
    hold on;
    str_legend = {};
    %% Different conventions
    for convention = 1:conventions_count
    plot(data_struct.MAP_D_grad_regular_interp(indices) * kBT,...
        (data_struct.MAP_fD_mean(lambda_type, indices, convention, 1) - data_struct.fD_theor_data(indices)'), markers_list{convention},...
        'markers', marker_size, 'LineWidth', line_width, 'color', color_sequence(convention, :));
        str_legend{end + 1} = conventions_names{convention};
    end;

    %% Identity-like lines (theory)
    % y = x
    h_theor_1 = plot(x_lim_vec, x_lim_vec, 'k--');
    % y = x/2
    h_theor_2 = plot(x_lim_vec, 1/2 * x_lim_vec, 'k--');
    % y = 0
    h_theor_0 = plot(x_lim_vec, 0 * x_lim_vec, 'k--');
    % y = -x/2
    h_theor_m2 = plot(x_lim_vec, -1/2 * x_lim_vec, 'k--');
    
    %% Adjust
    % Legend
    if lambda_type == 1
        legend(str_legend, 'location', 'north', 'FontSize', font_size-3);
    end;
    
    xlim(x_lim_vec);
    ylim(y_lim_vec);
    
    % Labels
    xlabel('$k_\mathrm{B}T\nabla D$', 'interpreter', 'latex');
    ylabel('$fD$ bias', 'interpreter', 'latex');
    title_str = {'$\lambda^* = 0$', '$\lambda^* = 0.5$', '$\lambda^* = 1$', 'Random $\lambda^*$'};
    title(title_str{lambda_type}, 'interpreter', 'latex');
    % Subplot label
    text(sublabel_x, sublabel_y, strcat('(', char('a' + lambda_type - 1), ')'), 'Units', 'Normalized', 'VerticalAlignment', 'Top');
    
    % Reorder curves
    uistack([h_theor_0, h_theor_m2, h_theor_2, h_theor_1], 'bottom');
end;


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










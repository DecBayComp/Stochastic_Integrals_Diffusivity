


function fig_count = plot_article_fD_bias(data_struct, fig_count, bl_save_figures)


%% Constants
sublabel_x = 0.02;
sublabel_y = 0.99;
load_constants;


%% === Plot fD bias as a function fo the gradient ===
fig_count = fig_count + 1; 
h_fig = figure(fig_count);
set_article_figure_size(h_fig, 1, 2, 1);
clf;
x_lim_vec = [-1, 1] * 0.17;
y_lim_vec = [-0.16, 0.175];
spacing = 0.08;
ML = 0.07;
MR = 0.02;
MT = 0.08;
MB = 0.15;
for l_ind = 1:lambda_count
    h_sub = subaxis(1, lambda_count, l_ind, 'Spacing', spacing, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);
    hold on;
    % Filter indices from one best period
    w = 10.0;
    x_left = (1/2 + 2*1)/w;
    x_right = (1/2 + 2*2)/w;
    indices = data_struct.x_bins_centers_saved{l_ind} >= x_left & data_struct.x_bins_centers_saved{l_ind} <= x_right;
    % Divine
    plot(data_struct.MAP_fwd_D_grad_regular_interp{l_ind}(indices) * kBT,...
        (data_struct.MAP_fwd_fD_divine_mean{l_ind}(1, indices) - data_struct.fD_theor_data{l_ind}(indices)'), markers_list{1}, 'markers', marker_size);
    % Ito
    plot(data_struct.MAP_fwd_D_grad_regular_interp{l_ind}(indices) * kBT,...
        (data_struct.MAP_fD_Ito_mean{l_ind}(1, indices) - data_struct.fD_theor_data{l_ind}(indices)'), markers_list{2}, 'markers', marker_size);
    % Stratonovich
    plot(data_struct.MAP_fwd_D_grad_regular_interp{l_ind}(indices) * kBT,...
        (data_struct.MAP_fwd_fD_Stratonovich_mean{l_ind}(1, indices) - data_struct.fD_theor_data{l_ind}(indices)'), markers_list{3}, 'markers', marker_size);
    % Marginalized
    plot(data_struct.MAP_fwd_D_grad_regular_interp{l_ind}(indices) * kBT,...
        (data_struct.MAP_fwd_fD_marginalized_mean{l_ind}(1, indices) - data_struct.fD_theor_data{l_ind}(indices)'), markers_list{4}, 'markers', marker_size);
    % Identity-like lines
    xlim(x_lim_vec);
    x_lim_vec = xlim();
    % y = x
    h_theor_1 = plot(x_lim_vec, x_lim_vec, 'k--');
        % y = x/2
    h_theor_2 = plot(x_lim_vec, 1/2 * x_lim_vec, 'k--');
    % y = 0
    h_theor_0 = plot(x_lim_vec, 0 * x_lim_vec, 'k--');
    % y = -x/2
    h_theor_m2 = plot(x_lim_vec, -1/2 * x_lim_vec, 'k--');
    % Legend
    if l_ind == 1
        str_legend = {'Orcl', 'Ito', 'Str', 'Mar'};
        legend(str_legend, 'location', 'southwest');
    end;
    % Add subplot label
    text(sublabel_x, sublabel_y, strcat('(', char('a' + l_ind - 1), ')'), 'Units', 'Normalized', 'VerticalAlignment', 'Top');
    % Reorder curves
    uistack([h_theor_0, h_theor_m2, h_theor_2, h_theor_1], 'bottom');
    xlabel('$k_\mathrm{B}T\nabla D$', 'interpreter', 'latex');
    ylabel('fD bias');
    title(sprintf('$\\lambda^* = %.2f$', lambda_array(l_ind)), 'interpreter', 'latex');
    % xlim([-0.5, 0.5]);
    ylim(y_lim_vec);
end;


%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize', [fig_pos(3), fig_pos(4)]);
% Set filename
output_filename = 'Force_bias_vs_gradient.pdf';
output_full_path = strcat(output_figures_folder, output_filename);
if bl_save_figures
    print(h_fig, output_full_path, '-dpdf', '-r0');
end;





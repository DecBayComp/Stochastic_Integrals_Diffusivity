

function fig_count = plot_article_D_posterior(data_struct, fig_count, bl_save_figures)



%% Constants
load_constants;
interval_scaling_factor_for_D = 6;
D_PRECISION = 1e-5;
D_steps = 1 + 2^6;
rows = 1;
cols = 3;


%% Plot diffusivity posterior in one bin, diffusivity profile and the fail rate
fig_count = fig_count + 1; 
h_fig = figure(fig_count);
set_article_figure_size(h_fig, rows, 2, 1);
clf;
% x_lim_vec = [D_min, D_max];
x_lim_vec = [7, 12] * 1e-3;
default_color_order = get(gca,'colororder');
% Initialize subplots
SH = 0.07;
SV = 0.125;
ML = 0.07;
MR = 0.02;
MT = 0.09;
MB = 0.18;
sublabel_x = 0.9;
sublabel_y = 0.1;
h_sub = subaxis(rows, cols, 1, 'SH', SH, 'SV', SV, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);
hold on;
% Get the selected bin
[selected_bins_indices, selected_bins_centers] = get_selected_bins_indices(data_struct);
%% Prepare mesh taking into account the most likely values
% Calculate the most probable values of D
D_values = zeros(1, lambda_count + 3);
for l_ind = 1:lambda_count 
    bin = selected_bins_indices(l_ind);
    [~, ~, nu_n, sigma2_n] = get_n_parameters(l_ind, bin, data_struct, 'forward');
    D_values(l_ind) = nu_n * sigma2_n / (2 * t_step * (2 + nu_n));
end;
% Calculate exact D
D_values(lambda_count + 1 : lambda_count + 3) = D_func(selected_D_case, selected_bins_centers * L, L);
% Detect the longest interval and centrally increase its length
D_min = min(D_values);
D_max = max(D_values);
D_interval = D_max - D_min;
D_min = max([D_min - D_interval * (interval_scaling_factor_for_D-1)/2, D_PRECISION]);    % D >= 0
D_max = D_max + D_interval * (interval_scaling_factor_for_D-1)/2;
% Prepare mesh
D_step = (D_max - D_min) / (D_steps - 1);
D_mesh = D_min:D_step:D_max;
% Add special points calculated earlier
D_mesh = sort([D_mesh, D_values]);

s1 = lambda_count;
s2 = length(D_mesh);
D_pdf_plot_data = zeros(s1, s2);
% D data
for l_ind = 1:lambda_count
    D_pdf_plot_data(l_ind, :) =  bin_D_posterior_func (l_ind,...
        selected_bins_indices(l_ind), D_mesh, t_step, data_struct, 'forward');
end;
% % Checking the PDF normalization
% pdf_norm = zeros(1, lambda_count);
% for lambda_ind = 1:lambda_count
%     pdf_norm(lambda_ind) = trapz(D_mesh, D_pdf_plot_data(lambda_ind, :));
% end;
% fprintf('Norm for D PDF:\n');
% disp(pdf_norm);
% Plot
for l_ind = 1:lambda_count
    plot(D_mesh, D_pdf_plot_data(l_ind, :), strcat('-', markers_list{l_ind}),...
        'color', default_color_order(l_ind, :), 'LineWidth', line_width);
end;
% Simulated value
y_lim_vec = [0, max(max(D_pdf_plot_data)) * 1.1];
plot(D_values(end) .* [1, 1], y_lim_vec, '--k', 'LineWidth', line_width);
ylim(y_lim_vec);    % I don't know why it changes scale without it
% Adjust
xlabel('D', 'interpreter', 'latex');
ylabel('PDF', 'interpreter', 'latex');
xlim(x_lim_vec);
title(sprintf('Posterior for $x^*\\approx%.2f$', selected_x_over_L), 'interpreter', 'latex');
% Add subplot label
text(sublabel_x, sublabel_y, '(a)', 'Units', 'Normalized', 'VerticalAlignment', 'Top');
box on;
% Legend
str_legend_local = {'$ \lambda^*=0$', '$ \lambda^*=1/2$', '$ \lambda^*=1$'};
legend(str_legend_local, 'location', 'northwest', 'interpreter', 'latex', 'fontsize', font_size);


%% === Plot diffusivity profile and fail rate ===
% fig_count = fig_count + 1; 
% h_fig = figure(fig_count);
% set_article_figure_size(h_fig, height_to_width_ratio_single_figure, 2);
% clf;
x_lim_vec = [0, x_max];

for l_ind = 1:lambda_count
        %% Profile plot
    subaxis(rows, cols, 2);
%     subplot(2, lambda_count, lambda_count + l_ind);
    hold on;
    % D
    plot(data_struct.x_bins_centers_saved{l_ind},  data_struct.MAP_D_mean{l_ind}(1, :),...
        strcat('-', markers_list{l_ind}), 'LineWidth', line_width);
%     data_struct.MAP_D_mean
    % Adjust
    xlim(x_lim_vec);
%     ylim([-0.1, 0.3]);
    ylim([0.004, 0.016]);
    xlabel('x', 'interpreter', 'latex');
    ylabel('D', 'interpreter', 'latex');
    title('Average D profile', 'interpreter', 'latex');
    % Add subplot label
    text(sublabel_x, sublabel_y, '(b)', 'Units', 'Normalized', 'VerticalAlignment', 'Top');
    box on;
    set(gca, 'LooseInset', get(gca,'TightInset'))
    
% % %     %% Success rate plot
% % %     h_sub = subaxis(rows, cols, 3);
% % %     hold on;
% % %     % D
% % %     plot(data_struct.x_bins_centers_saved{l_ind},  data_struct.UR_D{l_ind},...
% % %         'LineWidth', line_width);
% % %     % Plot the used confidence level
% % %     plot(x_lim_vec, [1, 1] * (1 - CONF_LEVEL), 'k--', 'linewidth', line_width);
% % %     % Adjust
% % % %     title(sprintf('$\\lambda^* = %.2f$', lambda_array(l_ind)), 'interpreter', 'latex');
% % %     xlim(x_lim_vec);
% % %     ylim([-0.02, 1.02]);
% % %     xlabel('x');
% % %     ylabel('Fail rate');
% % %     title('Fail rate');
% % %     box on;
% % %     pos = get(gca, 'Position');
% % %     set(gca, 'LooseInset', get(gca,'TightInset'))
    % Legend
%     if l_ind == 2
%         legend_position = [0.411, 0.583, 0.067, 0.127];
%         h_leg = legend(str_legend, 'position', legend_position);
%         legend boxon;
%     end;    
end;
subaxis(rows, cols, 2);
% Theory
plot(data_struct.x_fine_mesh, data_struct.D_theor_fine_data, '--k', 'LineWidth', line_width);






% %% === Average minimum confidence level including the true value ===
% y_lim_vec = [30, 101];
% for l_ind = 1:lambda_count
%     h_sub = subaxis(rows, cols, 4);
%     hold on;
%     % D
%     plot(data_struct.x_bins_centers_saved{l_ind},  data_struct.MAP_D_mean{l_ind}(4, :) * 100,...
%         'LineWidth', line_width);
%     % Adjust
% %     title(sprintf('$\\lambda^* = %.2f$', lambda_array(l_ind)), 'interpreter', 'latex');
%     xlim(x_lim_vec);
% %     y_lim_vec = ylim();
% %     y_lim_vec(2) = 100;
%     ylim(y_lim_vec);
%     xlabel('x');
%     ylabel('Av. min. CI');
%     title('Average minimum confidence interval');
%     box on;
% 
% end;


%% === Bias in D as a function of the D' within one period ===
x_lim_vec = [-1, 1] * 1.50e-3;
% y_lim_vec = [30, 101];
for l_ind = 1:lambda_count
    h_sub = subaxis(rows, cols, 3);
    hold on;
    % Filter indices from one best period
    w = 10.0;
    x_left = (1/2 + 2*1)/w;
    x_right = (1/2 + 2*2)/w;
    indices = data_struct.x_bins_centers_saved{l_ind} >= x_left & data_struct.x_bins_centers_saved{l_ind} <= x_right;
%     indices = data_struct.x_bins_centers_saved{l_ind} >= -1/2;
    % D
    difference = data_struct.MAP_D_mean{l_ind}(1, indices) - data_struct.D_theor_data{l_ind}(1, indices);
    half_error = (data_struct.MAP_D_mean{l_ind}(2, indices) + data_struct.MAP_D_mean{l_ind}(3, indices))/2;
    plot(difference, half_error, markers_list{l_ind}, 'markers', marker_size);
%     plot(data_struct.D_theor_data{l_ind}(3, indices),...
%         (data_struct.MAP_D_mean{l_ind}(1, indices) - data_struct.D_theor_data{l_ind}(indices)'), 'o');
% %     % It�
% %     plot(data_struct.MAP_fwd_D_grad_regular_interp{l_ind}(indices),...
% %         (data_struct.MAP_fD_Ito_mean{l_ind}(1, indices) - data_struct.fD_theor_data{l_ind}(indices)') / kBT, 'o', 'markers', marker_size);
% %     % Stratonovich
% %     plot(data_struct.MAP_fwd_D_grad_regular_interp{l_ind}(indices),...
% %         (data_struct.MAP_fwd_fD_Stratonovich_mean{l_ind}(1, indices) - data_struct.fD_theor_data{l_ind}(indices)') / kBT, '+', 'markers', marker_size);
% %     % Marginalized
% %     plot(data_struct.MAP_fwd_D_grad_regular_interp{l_ind}(indices),...
% %         (data_struct.MAP_fwd_fD_marginalized_mean{l_ind}(1, indices) - data_struct.fD_theor_data{l_ind}(indices)') / kBT, '^', 'markers', marker_size);
    % Identity
    y_lim_vec = ylim();
    x_temp = [ - y_lim_vec(2), 0, y_lim_vec(2)];
    y_temp = [-1, 0, 1] .* x_temp;
    plot(x_temp, y_temp, '--k');
end;
% Adjust
xlim(x_lim_vec);
% ylim(y_lim_vec);
xlabel('D bias', 'interpreter', 'latex');
ylabel('Half confidence interval', 'interpreter', 'latex');
title('Error cone', 'interpreter', 'latex');
% Add subplot label
text(sublabel_x, sublabel_y, '(c)', 'Units', 'Normalized', 'VerticalAlignment', 'Top');



%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize', [fig_pos(3), fig_pos(4)]);
% Set filename
output_filename = 'Diffusivity_profile_fail_rate.pdf';
output_full_path = strcat(output_figures_folder, output_filename);
if bl_save_figures
    print(h_fig, output_full_path, '-dpdf', '-r0');
end;


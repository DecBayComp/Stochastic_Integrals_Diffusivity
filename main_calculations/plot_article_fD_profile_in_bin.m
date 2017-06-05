


function fig_count = plot_article_fD_profile_in_bin(data_struct, fig_count, bl_save_figures)


%% Constants
load_constants;
fD_steps = 1 + 2^7;
x_lim_vec = [0, 0.41];


%% === Plot fD profile in one bin for each inference type ===
fig_count = fig_count + 1; 
h_fig = figure(fig_count);
set_article_figure_size(h_fig, 1, 1, 1);
clf;
hold on;
% l_ind = 3;
[selected_bins_indices, selected_bins_centers] = get_selected_bins_indices(data_struct);
% Calculate the most probable values of f*D
fD_values = zeros(1, lambda_count + 3);
for l_ind = 1:lambda_count 
    bin = selected_bins_indices(l_ind);
    [mu_n, ~, ~, ~] = get_n_parameters(l_ind, bin, data_struct, 'forward');
    fD_values(l_ind) = mu_n * kBT / t_step;
end;
%% Calculate exact fD
fD_values(lambda_count + 1: lambda_count + 3) = f_func(selected_f_case, selected_bins_centers*L, L)...
    .* D_func(selected_D_case, selected_bins_centers*L, L);
% Detect the longest interval and centrally increase its length
factor = 5;
fD_min = min(fD_values);
fD_max = max(fD_values);
fD_interval = fD_max - fD_min;
fD_min = fD_min - fD_interval * (factor-1)/2;
fD_max = fD_max + fD_interval * (factor-1)/2;
% Prepare mesh
fD_step = (fD_max - fD_min) / (fD_steps - 1);
fD_mesh = fD_min:fD_step:fD_max;
% Add special points
fD_mesh = sort([fD_mesh, fD_values]);
% Calculate data to plot
size1 = lambda_count;
size2 = length(fD_mesh);
fD_pdf_plot_data_divine = zeros(size1, size2);
fD_pdf_plot_data_Ito = zeros(size1, size2);
fD_pdf_plot_data_Stratonovich = zeros(size1, size2);
fD_pdf_plot_data_marginalized = zeros(size1, size2);
% fD data
for l_ind = 1:lambda_count
    bin = selected_bins_indices(l_ind);
    % Divine
    fD_pdf_plot_data_divine(l_ind, :) = bin_fD_divine_inference_posterior_func(data_struct, ...
        l_ind, bin, fD_mesh, data_struct.MAP_fwd_D_grad_regular_interp{l_ind}(bin), 'forward');
    % Ito
    fD_pdf_plot_data_Ito(l_ind, :) = bin_fD_posterior_func (data_struct, l_ind, bin, fD_mesh, 'forward');    
    % Stratonovich
    fD_pdf_plot_data_Stratonovich(l_ind, :) = bin_fD_simple_Stratonovich_posterior_func(data_struct, ...
                l_ind, bin, fD_mesh, data_struct.MAP_fwd_D_grad_regular_interp{l_ind}(bin), 'forward');
    % Marginalized
    fD_pdf_plot_data_marginalized(l_ind, :) = bin_fD_lambda_marginalized_posterior_func(data_struct, ...
                l_ind, bin, fD_mesh, data_struct.MAP_fwd_D_grad_regular_interp{l_ind}(bin), 'forward');
    
end;
% Plot
l_ind = 3;
% Divine
plot(fD_mesh, fD_pdf_plot_data_divine(l_ind, :), strcat('-', markers_list{1}), 'LineWidth', line_width/2);
% Ito
plot(fD_mesh, fD_pdf_plot_data_Ito(l_ind, :), strcat('-', markers_list{2}), 'LineWidth', line_width/2);
% Stratonovich
plot(fD_mesh, fD_pdf_plot_data_Stratonovich(l_ind, :), strcat('-', markers_list{3}), 'LineWidth', line_width/2);
% Marginalized
plot(fD_mesh, fD_pdf_plot_data_marginalized(l_ind, :), '-', 'LineWidth', line_width);
% True value
y_lim_vec = [0, max([fD_pdf_plot_data_marginalized(l_ind, :), fD_pdf_plot_data_Ito(l_ind, :)]) * 1.05];
plot(fD_values(lambda_count + 1) .* [1, 1], y_lim_vec, '--k', 'LineWidth', line_width);
% ylim(y_lim_vec);    % I don't know why it changes scale without it
% Adjust
title(sprintf('Posterior for $x\\approx%.2f$, $\\lambda^* = %.1f$', selected_x_over_L, lambda_array(l_ind)), 'interpreter', 'latex');
xlim(x_lim_vec);
ylim(y_lim_vec);
xlabel('fD', 'interpreter', 'latex');
ylabel('PDF', 'interpreter', 'latex');
% Legend
str_legend_local = {'Orcl', 'Ito', 'Str', 'Mar'};
legend(str_legend_local, 'location', 'northwest', 'interpreter', 'latex', 'fontsize', font_size);


%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize', [fig_pos(3), fig_pos(4)]);
% Set filename
output_filename = 'Force_posterior_in_one_bin.pdf';
output_full_path = strcat(output_figures_folder, output_filename);
if bl_save_figures
    print(h_fig, output_full_path, '-dpdf', '-r0');
end;
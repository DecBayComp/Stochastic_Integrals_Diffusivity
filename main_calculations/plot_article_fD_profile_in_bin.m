


function plot_article_fD_profile_in_bin(data_struct, trials_data, fig_count, bl_save_figures)
%% === Plot fD profile in one bin for each inference type ===


%% Constants
load_constants;
fD_steps = 1 + 2^7.5;
factor = 8;
y_factor = 1e3;
lambda_type = enum_lambda_rand;


%% Initialize
cur_data_struct = trials_data{data_struct.trial_first_simulation_type_index(lambda_type)};
lambda = cur_data_struct.lambda;
[selected_bins_indices, selected_bins_centers] = get_selected_bins_indices(cur_data_struct);
bin = selected_bins_indices;
D_grad = cur_data_struct.MAP_D_grad_regular_interp(bin);


%% Calculate
% The most probable values of fD for each convention
fD_values = zeros(1, conventions_count + 1);
[mu_n, ~, ~, ~] = get_n_parameters(bin, cur_data_struct, 'forward');
% Ito
tmp_lambda = 0;
fD_values(1) = (mu_n / t_step - tmp_lambda * D_grad) * kBT;
% Stratonovich
tmp_lambda = 0.5;
fD_values(2) = (mu_n / t_step - tmp_lambda * D_grad) * kBT;
% Hanggi
tmp_lambda = 1;
fD_values(3) = (mu_n / t_step - tmp_lambda * D_grad) * kBT;
% Marginalized
tmp_lambda = 0.5;
fD_values(4) = (mu_n / t_step - tmp_lambda * D_grad) * kBT;
% Oracle
tmp_lambda = lambda;
fD_values(5) = (mu_n / t_step - tmp_lambda * D_grad) * kBT;
% Exact fD
fD_values(conventions_count + 1) = f_func(selected_f_case, selected_bins_centers*L, L)...
    .* D_func(selected_D_case, selected_bins_centers*L, L);

% Among the calculated values, detect the longest interval and centrally increase its length
fD_min = min(fD_values);
fD_max = max(fD_values);
fD_interval = fD_max - fD_min;
fD_min = fD_min - fD_interval * (factor-1)/2;
fD_max = fD_max + fD_interval * (factor-1)/2;
% Prepare mesh
fD_step = (fD_max - fD_min) / (fD_steps - 1);
fD_mesh = fD_min:fD_step:fD_max;
% Add most likely values
fD_mesh = sort([fD_mesh, fD_values]);
fD_mesh_length = length(fD_mesh);


%% Calculate PDF values on mesh
fD_data = zeros(conventions_count, fD_mesh_length);
% Oracle
fD_data(enum_conv_divine, :) = bin_fD_divine_inference_posterior_func(cur_data_struct, ...
    lambda, bin, fD_mesh, D_grad, 'forward');
% Ito
fD_data(enum_conv_Ito, :) = bin_fD_posterior_func (cur_data_struct, bin, fD_mesh, 'forward');    
% Stratonovich
fD_data(enum_conv_Stratonovich, :) = bin_fD_simple_Stratonovich_posterior_func(cur_data_struct, ...
            bin, fD_mesh, D_grad, 'forward');
% Hanggi
fD_data(enum_conv_Hanggi, :) = bin_fD_simple_Hanggi_posterior_func(cur_data_struct, ...
    bin, fD_mesh, D_grad, 'forward');
% Marginalized
fD_data(enum_conv_marginalized, :) = bin_fD_lambda_marginalized_posterior_func(cur_data_struct, ...
            bin, fD_mesh, D_grad, 'forward');

        
%% Plot
h_fig = figure(fig_count);
set_article_figure_size(h_fig, 1, 1, 1);
clf;
hold on;
str_legend = {};
% PDFs
for convention = [enum_conv_Ito, enum_conv_Stratonovich, enum_conv_Hanggi, enum_conv_divine]
    plot(fD_mesh, fD_data(convention, :), strcat('-', markers_list{convention}), 'color', color_sequence(convention, :),...
        'LineWidth', line_width/2, 'markers', marker_size);
    str_legend{length(str_legend) + 1} = conventions_names{convention};
end
plot(fD_mesh, fD_data(enum_conv_marginalized, :), '-', 'color', color_sequence(enum_conv_marginalized, :),...
    'LineWidth', line_width*2, 'markers', marker_size);
str_legend{length(str_legend) + 1} = conventions_names{enum_conv_marginalized};

% Exact value
max_pdf = max(max(fD_data));
y_lim_vec = [0, max_pdf * 1.05];
h_theor = plot(fD_values(conventions_count + 1) .* [1, 1], y_lim_vec, '--k', 'LineWidth', line_width);

%% Adjust
% X limits
indices = sum(fD_data >= max_pdf / y_factor, 1);
fD_min = fD_mesh(find(indices, 1, 'first'));
fD_max = fD_mesh(find(indices, 1, 'last'));
x_lim_vec = [fD_min, fD_max];
xlim(x_lim_vec);
ylim(y_lim_vec);
box on;

xlabel('$f^AD$', 'interpreter', 'latex');
ylabel('PDF', 'interpreter', 'latex');
title(sprintf('$x\\approx%.2f$, $\\lambda^* = %.2f$', selected_x_over_L, lambda), 'interpreter', 'latex');
% Legend
legend(str_legend, 'location', 'northeast', 'interpreter', 'latex', 'fontsize', legend_font_size);
% Reorder curves
uistack(h_theor, 'bottom');


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








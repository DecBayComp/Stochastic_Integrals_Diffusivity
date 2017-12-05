


function plot_article_a_profile_in_bin(data_struct, trials_data, lambda_type, sel_x_over_L)
%% === Plot fD profile in one bin for each inference type ===


%% Constants
load_constants;
a_steps = round(1 + 2^8);
factor = 8;
y_factor = 1e3;
% output_filename = 'a_one_bin.pdf';


%% Initialize
cur_data_struct = trials_data{data_struct.trial_first_simulation_type_index(lambda_type)+4};
lambda = cur_data_struct.lambda;
[selected_bins_indices, selected_bins_centers] = get_selected_bins_indices(cur_data_struct, sel_x_over_L);
bin = selected_bins_indices;
bb_prime = cur_data_struct.MAP_bb_prime_regular_interp(bin);


%% Calculate
% The most probable values of fD for each convention
a_values = zeros(1, conventions_count + 1);
[mu_n, ~, ~, ~] = get_n_parameters(bin, cur_data_struct, 'forward');
% Ito
tmp_lambda = 0;
a_values(1) = mu_n / t_step - tmp_lambda * bb_prime;
% Stratonovich
tmp_lambda = 0.5;
a_values(2) = mu_n / t_step - tmp_lambda * bb_prime;
% Hanggi
tmp_lambda = 1;
a_values(3) = mu_n / t_step - tmp_lambda * bb_prime;
% Marginalized
tmp_lambda = 0.5;
a_values(4) = mu_n / t_step - tmp_lambda * bb_prime;
% Oracle
tmp_lambda = lambda;
a_values(5) = mu_n / t_step - tmp_lambda * bb_prime;
% Exact a
a_values(conventions_count + 1) = data_struct.a_theor_data(bin);

% Among the calculated values, detect the longest interval and centrally increase its length
a_min = min(a_values);
a_max = max(a_values);
a_interval = a_max - a_min;
a_min = a_min - a_interval * (factor-1)/2;
a_max = a_max + a_interval * (factor-1)/2;
% Prepare mesh
a_step = (a_max - a_min) / (a_steps - 1);
a_mesh = a_min:a_step:a_max;
% Add most likely values
a_mesh = sort([a_mesh, a_values]);
a_mesh_length = length(a_mesh);



%% Calculate PDF values on mesh
a_data = zeros(conventions_count, a_mesh_length);
% Oracle
a_data(enum_conv_divine, :) = bin_a_divine_inference_posterior_func(cur_data_struct, ...
    lambda, bin, a_mesh, bb_prime, 'forward');
% Ito
a_data(enum_conv_Ito, :) = bin_a_posterior_func (cur_data_struct, bin, a_mesh, 'forward');    
% Stratonovich
a_data(enum_conv_Stratonovich, :) = bin_a_simple_Stratonovich_posterior_func(cur_data_struct, ...
            bin, a_mesh, bb_prime, 'forward');
% Hanggi
a_data(enum_conv_Hanggi, :) = bin_a_simple_Hanggi_posterior_func(cur_data_struct, ...
    bin, a_mesh, bb_prime, 'forward');
% Marginalized
a_data(enum_conv_marginalized, :) = bin_a_lambda_marginalized_posterior_func(cur_data_struct, ...
            bin, a_mesh, bb_prime, 'forward');

		
        
%% Plot
% figure(h_fig);
% set_article_figure_size(h_fig, 1, 1, 1);
% clf;
hold on;
str_legend = {};
% PDFs
for convention = [enum_conv_Ito, enum_conv_Stratonovich, enum_conv_Hanggi, enum_conv_divine]
    plot(a_mesh, a_data(convention, :), strcat('-', markers_list{convention}), 'color', color_sequence(convention, :),...
        'LineWidth', line_width, 'markers', marker_size);
    str_legend{length(str_legend) + 1} = conventions_names{convention};
end
plot(a_mesh, a_data(enum_conv_marginalized, :), '-', 'color', color_sequence(enum_conv_marginalized, :),...
    'LineWidth', line_width+1, 'markers', marker_size);
str_legend{length(str_legend) + 1} = conventions_names{enum_conv_marginalized};

% Exact value
max_pdf = max(max(a_data));
y_lim_vec = [0, max_pdf * 1.05];
h_theor = plot(a_values(conventions_count + 1) .* [1, 1], y_lim_vec, '--k', 'LineWidth', line_width);

%% Adjust
% X limits
indices = sum(a_data >= max_pdf / y_factor, 1);
a_min = a_mesh(find(indices, 1, 'first'));
a_max = a_mesh(find(indices, 1, 'last'));
x_lim_vec = [a_min, a_max];
xlim(x_lim_vec);
ylim(y_lim_vec);
box on;

xlabel('$a$, $\mu \mathrm{m/s}$', 'interpreter', 'latex');
ylabel('PDF', 'interpreter', 'latex');
title(sprintf('$x\\approx%.2f\\ \\mu \\mathrm{m}$, $\\lambda^* = %.2f$', selected_bins_centers, lambda), 'interpreter', 'latex');
% Legend
legend(str_legend, 'location', 'northeast', 'interpreter', 'latex', 'fontsize', legend_font_size);
% Reorder curves
uistack(h_theor, 'bottom');


% % % %% Save figure
% % % % Prepare printer
% % % h_fig.PaperPositionMode = 'auto';
% % % h_fig.Units = 'Inches';
% % % fig_pos = h_fig.Position;
% % % set(h_fig, 'PaperUnits','Inches','PaperSize', [fig_pos(3), fig_pos(4)]);
% % % % Set filename
% % % output_full_path = strcat(output_figures_folder, output_filename);
% % % if bl_save_figures
% % %     print(h_fig, output_full_path, '-dpdf', '-r0');
% % % end;
% % % 
% % % selected_bins_centers






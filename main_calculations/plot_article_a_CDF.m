


function plot_article_a_CDF(data_struct, trials_data, lambda_type, sel_x_over_L)
%% === Plot fD profile in one bin for each inference type ===


%% Constants
load_constants;
A_steps = round(1 + 2^7);
factor = 8;
y_factor = 1e3;
output_filename = 'a_one_bin.pdf';


%% Initialize
cur_data_struct = trials_data{data_struct.trial_first_simulation_type_index(lambda_type)};
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
a_values(conventions_count + 1) = f_func(selected_f_case, selected_bins_centers*L, L) / gamma_drag;

% Among the calculated values, detect the longest interval and centrally increase its length
a_min = min(a_values);
a_max = max(a_values);
a_interval = a_max - a_min;
a_min = a_min - a_interval * (factor-1)/2;
a_max = a_max + a_interval * (factor-1)/2;

% In the interval take the max. abs value
A_max = max(abs([a_min, a_max]));

% Prepare mesh
A_step = A_max / (A_steps - 1);
A_mesh = 0:A_step:A_max;

% Add most likely values
A_mesh = sort([A_mesh, abs(a_values)]);
A_mesh_length = length(A_mesh);



%% Calculate CDF values on mesh
A_data = zeros(conventions_count, A_mesh_length);

% Oracle
pdf_wrap_func = @(a) bin_a_divine_inference_posterior_func(cur_data_struct, lambda, bin, a, bb_prime, 'forward');
cdf_wrap_func = @(A) integral(pdf_wrap_func, -A, A, 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE);
for i = 1:A_mesh_length
	A_data(enum_conv_divine, i) = 1 - cdf_wrap_func(A_mesh(i));
end;

% Ito
pdf_wrap_func = @(a) bin_a_posterior_func (cur_data_struct, bin, a, 'forward');  
cdf_wrap_func = @(A) integral(pdf_wrap_func, -A, A, 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE);
for i = 1:A_mesh_length
	A_data(enum_conv_Ito, i) = 1 - cdf_wrap_func(A_mesh(i));
end;

% Stratonovich
pdf_wrap_func = @(a) bin_a_simple_Stratonovich_posterior_func(cur_data_struct, bin, a, bb_prime, 'forward');
cdf_wrap_func = @(A) integral(pdf_wrap_func, -A, A, 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE);
for i = 1:A_mesh_length
	A_data(enum_conv_Stratonovich, i) = 1 - cdf_wrap_func(A_mesh(i));
end;		

% Hanggi
pdf_wrap_func = @(a) bin_a_simple_Hanggi_posterior_func(cur_data_struct, bin, a, bb_prime, 'forward');
cdf_wrap_func = @(A) integral(pdf_wrap_func, -A, A, 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE);
for i = 1:A_mesh_length
	A_data(enum_conv_Hanggi, i) = 1 - cdf_wrap_func(A_mesh(i));
end;

% Marginalized
pdf_wrap_func = @(a) bin_a_lambda_marginalized_posterior_func(cur_data_struct, bin, a, bb_prime, 'forward');
cdf_wrap_func = @(A) integral(pdf_wrap_func, -A, A, 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE);
A_tmp = zeros(A_mesh_length, 1);
parfor i = 1:A_mesh_length
	A_tmp(i) = 1 - cdf_wrap_func(A_mesh(i));
end;		
A_data(enum_conv_marginalized, :) = A_tmp';
		
        
%% Plot
% figure(h_fig);
% set_article_figure_size(h_fig, 1, 1, 1);
% clf;
hold on;
str_legend = {};
% PDFs
for convention = [enum_conv_Ito, enum_conv_Stratonovich, enum_conv_Hanggi, enum_conv_divine]
    plot(A_mesh, A_data(convention, :), '-', 'color', color_sequence(convention, :),...
        'LineWidth', line_width, 'markers', marker_size);
    str_legend{length(str_legend) + 1} = conventions_names{convention};
end
plot(A_mesh, A_data(enum_conv_marginalized, :), '-', 'color', color_sequence(enum_conv_marginalized, :),...
     'LineWidth', line_width+1, 'markers', marker_size);
str_legend{length(str_legend) + 1} = conventions_names{enum_conv_marginalized};

% Exact value
max_pdf = max(max(A_data));
y_lim_vec = [0, 1];
h_theor = plot(a_values(conventions_count + 1) .* [1, 1], y_lim_vec, '--k', 'LineWidth', line_width);

%% Adjust
% X limits
indices = sum(A_data >= max_pdf / y_factor, 1);
a_min = A_mesh(find(indices, 1, 'first'));
a_max = A_mesh(find(indices, 1, 'last'));
x_lim_vec = [a_min, a_max];
xlim(x_lim_vec);
ylim(y_lim_vec);
box on;

xlabel('$A$, $\mu \mathrm{m/s}$', 'interpreter', 'latex');
ylabel('$\mathrm{Pr}(|a|>A)$', 'interpreter', 'latex');
title(sprintf('$x\\approx%.2f\\ \\mu \\mathrm{m}$, $\\lambda^* = %.2f$', selected_bins_centers, lambda), 'interpreter', 'latex');
% Legend
legend(str_legend, 'location', 'northeast', 'interpreter', 'latex', 'fontsize', legend_font_size);
% Reorder curves
% uistack(h_theor, 'bottom');


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






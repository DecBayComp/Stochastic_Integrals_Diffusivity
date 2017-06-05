


function plot_article_inferred_D_profile(data_struct)
% This function performs inference of D in all bins for the selected f and
% D case


%% Globals


%% Constants
load_constants;
selected_lambda_case = 2;

height_to_width_ratio = 1/1.4;
output_figure_filename = 'Inferred_force_times_diffusivity_profile.pdf';
fD_PRECISION = 1e-7;
fD_ABS_MAX = 5; % 8 * 2 * 10;
D_GRAD_MAX = 30;
marker_size = 11;
fD_local_min = -0.5;
fD_local_max = 0.5;
% f_case_number = 1;
% D_case_number = 1;
% l_ind = 1;


%% Initialize
f_case_number = selected_f_case;
D_case_number = selected_D_case;
% Load data for the given D and f cases
filename = sprintf('D_%i_f_%i_data.mat', D_case_number, f_case_number);
full_path = strcat(output_data_folder, filename);
load(full_path, '-mat');

x_bin_widths = x_bins_widths_saved{selected_lambda_case};
small_shift = min([min(x_bin_widths) / 2, 0.015]) .* [1, -1];
% Legend
str_legend = cell(1, 1);
str_legend{1} = 'Simple It�';
str_legend{2} = 'Simple H�nggi';
str_legend{3} = 'Marginalized';
str_legend{4} = 'True';



%% Calculate D values and error bars
l_ind = selected_lambda_case;
x_bins_number = x_bins_number_saved(l_ind);
x_bins_centers = x_bins_centers_saved{l_ind};
inferred_MAP_D = zeros(3, x_bins_number);
for bin = 1:x_bins_number
    [~, ~, nu_n, sigma2_n] = get_n_parameters(l_ind, bin, data_struct, 'forward');
    inferred_MAP_D(1, bin) = nu_n / (nu_n + 2) * sigma2_n / (2 * t_step);
end;


%% Regularize D
[inferred_MAP_D_reg, inferred_MAP_D_grad_reg, inferred_MAP_D_grad_reg_interpolated, norm_cost, x_grad_mesh] =...
    regularize_gradient(inferred_MAP_D(1,:), x_bins_centers, alpha_reg);


%% Theoretical profile
x_min = min(x_bins_centers);
x_max = max(x_bins_centers);
x_theor_mesh = x_min:0.005:x_max;
[D_theor, D_grad_theor] = D_func(D_case_number, x_theor_mesh, L);


%% Plot
% Initialize figure
h_fig = figure(13);
set_article_figure_size(h_fig, height_to_width_ratio, 2);
default_color_order = get(gca,'colororder');
clf;
hold on;
% MLE
plot(x_bins_centers, inferred_MAP_D(1, :), 'b', 'LineWidth', 2);
% Regularized D
plot(x_bins_centers, inferred_MAP_D_reg, 'r', 'LineWidth', 2);
% Theory
plot(x_theor_mesh, D_theor, '--k', 'LineWidth', 1);

% Labels
xlabel('x');
ylabel('D');



%% === Calculate the gradient ===
%% Simple scheme with MAP values
x_bins_steps = x_bins_centers(2:end) - x_bins_centers(1:end - 1);
% x_grad_mesh = (x_bins_centers(2:end) + x_bins_centers(1:end - 1)) / 2;
simple_D_grad = (inferred_MAP_D(1, 2:end) - inferred_MAP_D(1, 1:end-1)) ./ x_bins_steps';


%% Calculate the real force posterior in one bin
fD_mesh = fD_local_min:0.01:fD_local_max;
bin = round(x_bins_number / 2);
% Marginalized force
fD_marginalized_posterior = bin_fD_lambda_marginalized_posterior_func(data_struct, l_ind, bin, fD_mesh, inferred_MAP_D_grad_reg_interpolated(bin), 'forward');
% Check norm
trapz(fD_mesh, fD_marginalized_posterior)
% A simple It� estimate
fD_simple_Ito = bin_fD_posterior_func (data_struct, l_ind, bin, fD_mesh, 'forward');
% % Simple Hanggi estimate
% tic;
% fprintf('Calculating the Hanggi estimate\n');
% fD_simple_Hanggi = bin_fD_simple_Hanggi_posterior_func(l_ind, bin, fD_mesh);
% fprintf('Finished: Calculating the Hanggi estimate. Elapsed: %.2fs\n', toc);
% Divine regularized inference
fD_divine = bin_fD_divine_inference_posterior_func(data_struct, l_ind, bin, fD_mesh, inferred_MAP_D_grad_reg_interpolated(bin), 'forward');
% Check norm
trapz(fD_mesh, fD_divine)




%% === Plot the gradeint ===
h_fig_grad = figure(14);
set_article_figure_size(h_fig_grad, height_to_width_ratio, 2);
clf;
hold on;
% Plot theory
plot(x_theor_mesh, D_grad_theor, 'k--', 'LIneWidth', 2);
% Simple difference
plot(x_grad_mesh, simple_D_grad, 'bo', 'LineWidth', 2)
% Regularized gradient
plot(x_grad_mesh, inferred_MAP_D_grad_reg, 'rx', 'LineWidth', 2);
% Regularized interpolated gradient
plot(x_bins_centers, inferred_MAP_D_grad_reg_interpolated, 'g-', 'LineWidth', 2);
% Adjust plot
xlabel('x');
ylabel('D''');
y_lim_vec = ylim();
y_lim_vec(1) = -0.75;
% ylim(y_lim_vec);


%% === Plot force posterior in one bin ===
h_fig_grad = figure(15);
set_article_figure_size(h_fig_grad, height_to_width_ratio, 2);
clf;
hold on;
% Marginalized force
plot(fD_mesh, fD_marginalized_posterior, 'LineWidth', 2);
% Simple It� inference
plot(fD_mesh, fD_simple_Ito, 'r-', 'LineWidth', 2);
% % Simple Hanggi inference
% plot(fD_mesh, fD_simple_Hanggi, 'g-', 'LineWidth', 2);
% Divine regularized inference
plot(fD_mesh, fD_divine, 'm-', 'LineWidth', 2);
% Real value
% % max_y = max(
fD_theor_value = f_func(f_case_number, x_bins_centers(bin)) * D_func(D_case_number, x_bins_centers(bin));
y_mesh = [0, max(fD_divine) * 1.05];
plot([1, 1] * fD_theor_value, y_mesh, 'k--');
% Adjust plot
xlabel('f*D');
ylabel('PDF');
% Legend
str_legend = cell(1, 3);
str_legend{1} = 'Marginalized force';
str_legend{2} = 'Simple It� inference';
% str_legend{3} = 'Simple Hanggi inference';
str_legend{3} = 'Regularized divine inference';
legend(str_legend, 'location', 'best');

















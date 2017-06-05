


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
str_legend{1} = 'Simple Itô';
str_legend{2} = 'Simple Hänggi';
str_legend{3} = 'Marginalized';
str_legend{4} = 'True';
% Initialize figure
h_fig = figure(13);
set_article_figure_size(h_fig, height_to_width_ratio, 2);
default_color_order = get(gca,'colororder');
clf;
hold on;


%% Theoretical equilibrium distribution
% x_theor_mesh = x_min:0.005:x_max;
% fD_data = D_func(D_case_number, x_theor_mesh, L) .* f_func(f_case_number, x_theor_mesh, L);


%% Calculate values and error bars
l_ind = selected_lambda_case;
x_bins_number = x_bins_number_saved(l_ind);
x_bins_centers = x_bins_centers_saved{l_ind};
inferred_MAP_fD_Ito = zeros(3, x_bins_number);
inferred_MAP_fD_Hanggi = zeros(3, x_bins_number);
inferred_MAP_fD_marginalized = zeros(3, x_bins_number);
% Choose bin
for bin = 2:x_bins_number
    fprintf('Calculating bin: %i/%i\n', bin, x_bins_number);
    
    % Calculate INFERENCE WITHOUT MARGINALIZATION, LAMBDA = 0
    fprintf('Calculating ITO inference\n');
    function_to_minimze = @(fD) bin_fD_posterior_func (l_ind, bin, fD, t_step, kBT);
    fD_Ito_inference = find_confidence_interval(function_to_minimze, [- fD_ABS_MAX, fD_ABS_MAX], true);
    % Save
    inferred_MAP_fD_Ito(:, bin) = fD_Ito_inference;
    
    % Calculate INFERENCE WITHOUT MARGINALIZATION, LAMBDA = 1
    fprintf('Calculating HÄNGGI inference\n');
    % Make Ito inference of fD in this bin
    function_to_minimze = @(fD) bin_fD_simple_Hanggi_posterior_func(l_ind, bin, fD);
%         optim_options = optimset('Display','iter', 'TolX', D_PRECISION);
%     optim_options = optimset('TolX', fD_PRECISION);
%     fD_MLE_Hanggi_inference = fminbnd(function_to_minimze, - fD_ABS_MAX, fD_ABS_MAX, optim_options);
    fD_MLE_Hanggi_inference = find_confidence_interval(function_to_minimze, [- fD_ABS_MAX, fD_ABS_MAX], false);
    % Save
    inferred_MAP_fD_Hanggi(:, bin) = fD_MLE_Hanggi_inference;
    
    % Calculate MARGINALIZED INFERENCE
    fprintf('Calculating MARGINALIZED inference\n');
    function_to_minimze = @(fD) bin_fD_marginalized_posterior_func (l_ind, bin, fD);
%     fD_MLE_marginalized = fminbnd(function_to_minimze, - fD_ABS_MAX, fD_ABS_MAX, optim_options);
    fD_MLE_marginalized = find_confidence_interval(function_to_minimze, [- fD_ABS_MAX, fD_ABS_MAX], false);
    % Save
    inferred_MAP_fD_marginalized(:, bin) = fD_MLE_marginalized;
    

    %% Plot
    clf;
    hold on;
    % Ito inference
    errorbar(x_bins_centers(2:bin), inferred_MAP_fD_Ito(1, 2:bin),...
            inferred_MAP_fD_Ito(2, 2:bin), inferred_MAP_fD_Ito(3, 2:bin), '-o', 'LineWidth', 1.5, 'color', 'b', 'markers', marker_size);
    % Hanggi inference
    errorbar(x_bins_centers(2:bin) + small_shift(1), inferred_MAP_fD_Hanggi(1, 2:bin),...
            inferred_MAP_fD_Hanggi(2, 2:bin), inferred_MAP_fD_Hanggi(3, 2:bin), '-o', 'LineWidth', 1.5, 'color', 'g', 'markers', marker_size);
    % Marginalized inference
    errorbar(x_bins_centers(2:bin) + small_shift(2), inferred_MAP_fD_marginalized(1, 2:bin),...
            inferred_MAP_fD_marginalized(2, 2:bin), inferred_MAP_fD_marginalized(3, 2:bin), '-o', 'LineWidth', 1.5, 'color', 'r', 'markers', marker_size);
    % Theory
    x_theor_min = x_bins_centers(2);
    x_theor_max = x_bins_centers(end);
    theor_steps_count = 1000;
    x_theor_step = (x_theor_max - x_theor_min) / (theor_steps_count - 1);
    x_theor_mesh = x_theor_min:x_theor_step:x_theor_max;
    fD_theor_data = D_func(D_case_number, x_theor_mesh, L) .* f_func(f_case_number, x_theor_mesh, L);
    plot(x_theor_mesh, fD_theor_data, '--k', 'LineWidth', 2);
    %     ylim([-30,30]);
    legend(str_legend, 'location', 'best');
    pause(0.01);
end;


%% Extract MLE diffusivity values



%% Plot
clf;
hold on;
% Ito inference
errorbar(x_bins_centers(2:x_bins_number), inferred_MAP_fD_Ito(1, 2:x_bins_number),...
        inferred_MAP_fD_Ito(2, 2:x_bins_number), inferred_MAP_fD_Ito(3, 2:x_bins_number), '-o', 'LineWidth', 1, 'color', 'b', 'markers', marker_size);
% % Hanggi inference
% errorbar(x_bins_centers(2:x_bins_number)+ small_shift(1), inferred_MAP_fD_Hanggi(1, 2:x_bins_number),...
%         inferred_MAP_fD_Hanggi(2, 2:x_bins_number), inferred_MAP_fD_Hanggi(3, 2:x_bins_number), '-o', 'LineWidth', 1, 'color', 'g', 'markers', marker_size);
% Marginalized inference
errorbar(x_bins_centers(2:x_bins_number) + small_shift(2), inferred_MAP_fD_marginalized(1, 2:x_bins_number),...
        inferred_MAP_fD_marginalized(2, 2:x_bins_number), inferred_MAP_fD_marginalized(3, 2:x_bins_number), '-o', 'LineWidth', 1, 'color', 'r', 'markers', marker_size);
% Theory
x_theor_min = x_bins_centers(2);
x_theor_max = x_bins_centers(end);
theor_steps_count = 1000;
x_theor_step = (x_theor_max - x_theor_min) / (theor_steps_count - 1);
x_theor_mesh = x_theor_min:x_theor_step:x_theor_max;
fD_theor_data = D_func(D_case_number, x_theor_mesh, L) .* f_func(f_case_number, x_theor_mesh, L);
plot(x_theor_mesh, fD_theor_data, '--k', 'LineWidth', 2);
%% Adust plot
xlabel('x');
ylabel('f*D');
box on;
% xlim([-0.5, 0.5]);
%         ylim([-fD_ABS_MAX, fD_ABS_MAX]);
ylim([- 0.1, 0.3]);
hold off;
legend(str_legend, 'location', 'best');


%% Save calculated data and intervals
% filename = sprintf('D_%i_f_%i_Profiles_with_confidence_intervals_data.mat', D_case_number, f_case_number);
% output_full_path = strcat(output_data_folder, filename);
% fprintf('Saving confidence intervals...\n');
% save(output_full_path, 'x_bins_centers', 'inferred_MAP_fD_simple', 'inferred_MAP_fD_marginalized', 'x_theor_mesh', 'fD_data');
% fprintf('Confidence intervals saved successfully!\n');


%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize',[fig_pos(3), fig_pos(4)]);

% Set filename
output_full_path = strcat(output_figures_folder, output_figure_filename);
print(h_fig, output_full_path, '-dpdf', '-r0');















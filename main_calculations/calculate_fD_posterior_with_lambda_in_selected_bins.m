

%% Globals
global selected_x_over_L;
global fD_marginalized_steps;
% global lambda_count;

%% Constants
interval_scaling_factor = 15;


%% Initialization
[selected_bins_indices, selected_bins_centers] = get_selected_bins_indices();

figure(8);
clf;
color_codes = get(gca,'colororder');

fprintf('Calculating marginalized force posterior...\n');



%% Calculate the most probable values of marginalized f*D to construct the mesh
fD_values = zeros(1, 3 * lambda_count + 3);
for lambda_ind = 1:lambda_count 
    bin = selected_bins_indices(lambda_ind);
    if bin ==1
        error('Cannot perform calculations in bin 1');
    end;
    x_step = x_bins_centers_saved{lambda_ind}(bin) - x_bins_centers_saved{lambda_ind}(bin - 1);
    % Bin j-1
    [mu_n1, kappa_n1, nu_n1, sigma2_n1] = get_n_parameters(lambda_ind, bin - 1);
    % Bin j
    [mu_n2, kappa_n2, nu_n2, sigma2_n2] = get_n_parameters(lambda_ind, bin);
    for i =1:3
        cur_lambda = (i-1) * 1/2;
        fD_values((lambda_ind - 1) * 3 + i) = kBT / t_step * (mu_n2 - cur_lambda / 2 / x_step...
            * (nu_n2 * sigma2_n2 / (3 + nu_n2) - nu_n1 * sigma2_n1 / (2 + nu_n1)));
    end;
end;

%% Calculate exact fD
fD_values(3 * lambda_count + 1 : 3 * lambda_count + 3) =...
    f_func(f_case_number, selected_bins_centers*L, L)...
    .* D_func(D_case_number, selected_bins_centers*L, L);

%% Detect the longest interval and centrally increase its length
fD_min = min(fD_values);
fD_max = max(fD_values);
fD_interval = fD_max - fD_min;
fD_min = fD_min - fD_interval * (interval_scaling_factor-1)/2;
fD_max = fD_max + fD_interval * (interval_scaling_factor-1)/2;

%% Prepare mesh
fD_step = (fD_max - fD_min) / (fD_marginalized_steps - 1);
fD_mesh = fD_min:fD_step:fD_max;
% Add special points
fD_mesh = sort([fD_mesh, fD_values]);

size1 = lambda_count;
fD_mesh_length = length(fD_mesh);
fD_pdf_plot_data = zeros(size1, fD_mesh_length);

selected_bins_count = length(selected_x_over_L);



% Evaluate integral on the mesh
    for lambda_ind = 1:lambda_count
        bin = selected_bins_indices(lambda_ind, 1);
        for fD_ind = 1:fD_mesh_length
            fD_pdf_plot_data(lambda_ind, fD_ind) =...
                bin_fD_marginalized_posterior_func (lambda_ind, bin, fD_mesh(fD_ind));



%% Plotting
fig_hand = figure(8);
set_my_fig_size (fig_hand);
clf;
hold on;
legends_str = cell(1, selected_bins_count);
% Solid lines (Ito)
for i_bin2 = 1:selected_bins_count
    lambda_ind2 = 1;
    plot(fD_mesh, squeeze(fD_pdf_plot_data(lambda_ind2, :)), '-',...
        'Color', color_codes(1, :), 'LineWidth', 2);
    legends_str = sprintf('x/L = %.3f', selected_bins_centers(lambda_ind2, 1)/L);
end;
legend(legends_str);
% Other lines
for i_bin2 = 1:selected_bins_count
    lambda_ind2 = 2;
    plot(fD_mesh, squeeze(fD_pdf_plot_data(lambda_ind2, :)), '--',...
        'Color', color_codes(1, :), 'LineWidth', 2);
    lambda_ind2 = 3;
    plot(fD_mesh, squeeze(fD_pdf_plot_data(lambda_ind2, :)), ':',...
        'Color', color_codes(1, :), 'LineWidth', 3);
end;


xlabel('f*D');
ylabel('PDF');
title('Marginalized posterior on f*D. \lambda=0 (solid), \lambda=0.5 (dashed), \lambda=1 (dotted)');
xlim([fD_min, fD_max]);
hold off;
pause(0.01);


        end;
        pause(0.01);
    end;


fprintf('Calculations completed!');
    
%% Checking the PDF normalization
pdf_norm = zeros(1, lambda_count);
for lambda_ind = 1:lambda_count
    pdf_norm(lambda_ind) = trapz(fD_mesh, fD_pdf_plot_data(lambda_ind, :));
end;
pdf_norm

%% Saving
output_filename_no_ext = sprintf('D%i_F%i_08_Posterior on force times diffusivity marginalized against lambda.png', D_case_number, f_case_number);
output_full_path_no_ext = strcat(output_figures_folder, output_filename_no_ext);
% output_full_path_pdf = strcat(output_figures_folder, output_filename_pdf);
if bl_save_figures
    fig_hand.PaperPositionMode = 'auto';
    saveas(fig_hand, output_full_path_no_ext, 'png');
    savefig(fig_hand, strcat(output_full_path_no_ext, '.fig'));
%     saveas(fig_hand, output_full_path_pdf, 'pdf');
    disp('Figure saved (force times diffusivity marginalized against lambda)');
end;










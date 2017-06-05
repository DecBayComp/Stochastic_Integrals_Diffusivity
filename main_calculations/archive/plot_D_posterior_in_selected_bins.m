

%% Globals
% global selected_x_over_L_coordinates;


%% Constants
D_steps = 1 + 2^10;


%% Choosing limits independently for each of the regimes
if D_case_number == 2 && f_case_number == 5
    D_min = 1e-3;
    D_max = 2;    
elseif D_case_number == 2 
    D_min = 1e-3;
    D_max = 1.5;   
elseif D_case_number == 3 && f_case_number == 2
    D_min = 1e-3;
    D_max = 0.1;  
elseif D_case_number == 3 && f_case_number == 3
    D_min = 1e-3;
    D_max = 0.12;
elseif D_case_number == 3 && f_case_number == 4
    D_min = 1e-3;
    D_max = 0.03;
elseif D_case_number == 3 && f_case_number == 5
    D_min = 1e-3;
    D_max = 0.2;
elseif D_case_number == 3
    D_min = 1e-3;
    D_max = 1;    
else
    D_min = 1e-3;
    D_max = 0.8;
end;


%% Initialization
[selected_bins_indices, selected_bins_centers] = get_selected_bins_indices();
selected_bins_count = size(selected_bins_centers, 2);


%% Preparing the mesh and data for plotting
% Mesh
D_step = (D_max - D_min) / (D_steps - 1);
D_mesh = D_min:D_step:D_max;
s1 = lambda_count;
s2 = D_steps;
D_pdf_plot_data = zeros(selected_bins_count, s1, s2);
% lambda_mesh = (1:lambda_count); %' * ones(1, s2);
% % % b_mesh = ones(s1, 1) * b_mesh;
% Evaluating
for i_bin = 1:selected_bins_count
    for lambda_ind = 1:lambda_count
        D_pdf_plot_data(i_bin, lambda_ind, :) = bin_D_pdf_func (lambda_ind,...
            selected_bins_indices(lambda_ind, i_bin),...
            D_mesh);
    end;
end;


%% Plotting
fig_hand = figure(5);
set_my_fig_size (fig_hand);
clf;
hold all;
color_codes = get(gca,'colororder');
legends_str = cell(1, selected_bins_count);
% Solid lines (Ito)
for i_bin = 1:selected_bins_count
    lambda_ind = 1;
    plot(D_mesh, squeeze(D_pdf_plot_data(i_bin, lambda_ind, :)), '-',...
        'Color', color_codes(i_bin, :), 'LineWidth', 2);
    legends_str{i_bin} = sprintf('x/L = %.3f', selected_bins_centers(l_ind, i_bin)/L);
end;
legend(legends_str);
% Other lines
for i_bin = 1:selected_bins_count
    lambda_ind = 2;
    plot(D_mesh, squeeze(D_pdf_plot_data(i_bin, lambda_ind, :)), '--',...
        'Color', color_codes(i_bin, :), 'LineWidth', 2);
    lambda_ind = 3;
    plot(D_mesh, squeeze(D_pdf_plot_data(i_bin, lambda_ind, :)), ':',...
        'Color', color_codes(i_bin, :), 'LineWidth', 3);
end;

xlabel('D');
ylabel('PDF');
title('Posterior on diffusivity in selected bins. \lambda=0 (solid), \lambda=0.5 (dashed), \lambda=1 (dotted)');
xlim([D_min, D_max]);
hold off;

%% Saving
output_filename_png = sprintf('D%i_F%i_05_Posterior on diffusivity.png', D_case_number, f_case_number);
output_full_path_svg = strcat(output_figures_folder, output_filename_png);
% output_full_path_pdf = strcat(output_figures_folder, output_filename_pdf);
if bl_save_figures
    fig_hand.PaperPositionMode = 'auto';
    saveas(fig_hand, output_full_path_svg, 'png');
%     saveas(fig_hand, output_full_path_pdf, 'pdf');
    disp('Figure saved (diffusivity posterior)');
end;











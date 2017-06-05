

%% Globals
% global selected_x_over_L_coordinates;

%% Constants
D_grad_steps = 1 + 2^9;

%% Choosing limits independently for each of the regimes
if D_case_number == 1 && f_case_number == 1
    D_grad_min = -15;
    D_grad_max = 10;
elseif D_case_number == 1 && f_case_number == 2
    D_grad_min = -0.4;
    D_grad_max = 0.75;   
elseif D_case_number == 1 && f_case_number == 3
    D_grad_min = -4;
    D_grad_max = 7;        
elseif D_case_number == 2 && f_case_number == 1
    D_grad_min = -40;
    D_grad_max = 20;
elseif D_case_number == 2 && f_case_number == 2
    D_grad_min = -3;
    D_grad_max = 6;    
elseif D_case_number == 2 && f_case_number == 3
    D_grad_min = -8;
    D_grad_max = 20;     
elseif D_case_number == 2 && f_case_number == 4
    D_grad_min = -30;
    D_grad_max = 30; 
elseif D_case_number == 2 && f_case_number == 5
    D_grad_min = -20;
    D_grad_max = 30;   
elseif D_case_number == 3 && f_case_number == 1
    D_grad_min = -0.05;
    D_grad_max = 0.2;
elseif D_case_number == 3 && f_case_number == 2
    D_grad_min = -0.03;
    D_grad_max = 0.12;
elseif D_case_number == 3 && f_case_number == 3
    D_grad_min = -0.05;
    D_grad_max = 0.15;
elseif D_case_number == 3 && f_case_number == 4
    D_grad_min = -0.01;
    D_grad_max = 0.06;
elseif D_case_number == 3 && f_case_number == 5
    D_grad_min = 0;
    D_grad_max = 0.01;
elseif D_case_number == 3
    D_grad_min = -10;
    D_grad_max = 10;    
else
    val = 8;
    D_grad_min = -val;
    D_grad_max = val;
end;




%% Initialization
[selected_bins_indices, selected_bins_centers] = get_selected_bins_indices();


%% Preparing the mesh and data for plotting
% Mesh
D_grad_step = (D_grad_max - D_grad_min) / (D_grad_steps - 1);
D_grad_mesh = D_grad_min:D_grad_step:D_grad_max;
size1 = lambda_count;
size2 = D_grad_steps;
D_grad_pdf_plot_data = zeros(selected_bins_count, size1, size2);
% Evaluating
for i_bin = 1:selected_bins_count
    for lambda_ind = 1:lambda_count
        D_grad_pdf_plot_data(i_bin, lambda_ind, :) = ...
            bin_D_grad_pdf_func (lambda_ind, selected_bins_indices(lambda_ind, i_bin),...
            D_grad_mesh);
    end;
end;


%% Plotting
fig_hand = figure(6);
set_my_fig_size (fig_hand);
clf;
hold all;
color_codes = get(gca,'colororder');
legends_str = cell(1, selected_bins_count);
% Solid lines (Ito)
for i_bin = 1:selected_bins_count
    lambda_ind = 1;
    plot(D_grad_mesh, squeeze(D_grad_pdf_plot_data(i_bin, lambda_ind, :)), '-',...
        'Color', color_codes(i_bin, :), 'LineWidth', 2);
    legends_str{i_bin} = sprintf('x/L = %.3f', selected_bins_centers(lambda_ind, i_bin)/L);
end;
legend(legends_str);
% Other lines
for i_bin = 1:selected_bins_count
    lambda_ind = 2;
    plot(D_grad_mesh, squeeze(D_grad_pdf_plot_data(i_bin, lambda_ind, :)), '--',...
        'Color', color_codes(i_bin, :), 'LineWidth', 2);
    lambda_ind = 3;
    plot(D_grad_mesh, squeeze(D_grad_pdf_plot_data(i_bin, lambda_ind, :)), ':',...
        'Color', color_codes(i_bin, :), 'LineWidth', 3);
end;

xlabel('grad D');
ylabel('PDF');
title('Posterior on the diffusivity gradient in selected bins. \lambda=0 (solid), \lambda=0.5 (dashed), \lambda=1 (dotted)');
xlim([D_grad_min, D_grad_max]);
hold off;

%% Saving
output_filename_png = sprintf('D%i_F%i_06_Posterior on diffusivity gradient.png', D_case_number, f_case_number);
output_full_path_svg = strcat(output_figures_folder, output_filename_png);
% output_full_path_pdf = strcat(output_figures_folder, output_filename_pdf);
if bl_save_figures
    fig_hand.PaperPositionMode = 'auto';
    saveas(fig_hand, output_full_path_svg, 'png');
%     saveas(fig_hand, output_full_path_pdf, 'pdf');
    disp('Figure saved (diffusivity gradient posterior)');
end;











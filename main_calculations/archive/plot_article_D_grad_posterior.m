



%% Globals
% global n_j;
% global V_j;

%% Constants
load_constants;
height_to_width_ratio = 1;
marker_step = 10;
D_grad_steps = 1 + 2^5;     % 2^9
D_grad_min = -2;
D_grad_max = 2;
interval_factor = 30;



%% Initialize figure
h_fig = figure(7);
clf;
set_article_figure_size(h_fig, height_to_width_ratio, 2);
default_color_order = get(gca,'colororder');
markers_list = {'-o','-s','-d','-^','-v'};


count = 1;
for f_case_number = 1:max_f_case_number
    for D_case_number = 1:max_D_case_number

        %% Load data
        filename = sprintf('D_%i_f_%i_data.mat', D_case_number, f_case_number);
        full_path = strcat(output_data_folder, filename);
        load(full_path, '-mat');
        
        %% Initialize
        subplot(max_f_case_number, max_D_case_number, count);
        hold on;
        [selected_bins_indices, selected_bins_centers] = get_selected_bins_indices();
        selected_bins_count = 1;
        
        %% Calculate the most probable values of D'
        D_grad_values = zeros(1, lambda_count + 3);
        for lambda_ind = 1:lambda_count 
            bin = selected_bins_indices(lambda_ind);
            x_step = x_bins_centers_saved{lambda_ind}(bin) - x_bins_centers_saved{lambda_ind}(bin - 1);
            [mu_n1, kappa_n1, nu_n1, sigma2_n1] = get_n_parameters(lambda_ind, bin - 1);
            [mu_n2, kappa_n2, nu_n2, sigma2_n2] = get_n_parameters(lambda_ind, bin);
            D_grad_values(lambda_ind) = (nu_n2 * sigma2_n2 /  (2 + nu_n2)...
                - nu_n1 * sigma2_n1 / (2 + nu_n1)) / (2 * t_step * x_step);
        end;
        
        %% Calculate exact D'
        [~, D_grad_val] = D_func(D_case_number, selected_bins_centers*L, L);
        D_grad_values(lambda_count+1 : lambda_count + 3) = D_grad_val;
        
        %% Detect the longest interval and centrally increase its length
        D_grad_min = min(D_grad_values);
        D_grad_max = max(D_grad_values);
        D_grad_interval = D_grad_max - D_grad_min;
%         if f_case_number ==3 && D_case_number == 3    % Manual override since the gradient is probably 0
%             val = 8;
%             D_grad_min = -val;
%             D_grad_max = val;
%         else
            D_grad_min = D_grad_min - D_grad_interval * (interval_factor-1)/2;
            D_grad_max = D_grad_max + D_grad_interval * (interval_factor-1)/2;
%         end;
        
        %% Prepare mesh
        D_grad_step = (D_grad_max - D_grad_min) / (D_grad_steps - 1);
        D_grad_mesh = D_grad_min:D_grad_step:D_grad_max;
        % Add special points
        D_grad_mesh = sort([D_grad_mesh, D_grad_values]);
        
        
        size1 = lambda_count;
        size2 = length(D_grad_mesh);
        D_grad_pdf_plot_data = zeros(size1, size2);
        % fD data
        for lambda_ind = 1:lambda_count
            bin = selected_bins_indices(lambda_ind);
            D_grad_pdf_plot_data(lambda_ind, :) = bin_D_grad_posterior_func (lambda_ind, bin, D_grad_mesh,...
                t_step);
        end;
        
        
        %% Checking the PDF normalization
        pdf_norm = zeros(1, lambda_count);
        for lambda_ind = 1:lambda_count
            pdf_norm(lambda_ind) = trapz(D_grad_mesh, D_grad_pdf_plot_data(lambda_ind, :));
        end;
        fprintf('Norm for F %i/%i, D %i/%i:\n', f_case_number, max_f_case_number, D_case_number, max_D_case_number);
        disp(pdf_norm);


        %% Plot
        for lambda_ind = 1:lambda_count
            plot(D_grad_mesh, D_grad_pdf_plot_data(lambda_ind, :), '-',...
                'color', default_color_order(lambda_ind, :), 'LineWidth', 2);
        end;


        %% Simulated value
        y_lim_vec = [0, max(max(D_grad_pdf_plot_data)) * 1.1];
        plot(D_grad_val .* [1, 1], y_lim_vec, '--k', 'LineWidth', 2);
        ylim(y_lim_vec);    % I don't know why it changes scale without it
% 
% 
        %% Adjust
        if f_case_number == 1
            title(sprintf('D%i', D_case_number), 'FontWeight','Bold');
        end;
        if f_case_number == max_f_case_number
            xlabel('$D_j''$', 'interpreter', 'latex');
        end;
        if D_case_number == 1
            ylabel(sprintf('F%i', f_case_number), 'FontWeight','Bold');
        end;
%         ylabel('PDF');
%         x_lim_vec = xlim();
%         x_lim_vec(2) = x_lim_vec(2) * 1.1;
%         xlim(x_lim_vec);
        x_lim_vec = [D_grad_min, D_grad_max];
        xlim(x_lim_vec);
        box on;
        
        count = count + 1;
        pause(0.01);
    end;
end;




%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize',[fig_pos(3), fig_pos(4)]);

% Set filename
output_filename = 'Disffusivity_gradient_posterior.pdf';
output_full_path = strcat(output_figures_folder, output_filename);
print(h_fig, output_full_path, '-dpdf', '-r0');








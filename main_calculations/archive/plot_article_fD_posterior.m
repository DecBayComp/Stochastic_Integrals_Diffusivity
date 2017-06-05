



%% Globals
load_constants;
% global n_j;
% global V_j;

%% Constants
height_to_width_ratio = 1;
f_case_number = 1;
D_case_number = 1;
marker_step = 10;
fD_steps = 1 + 2^10;
fD_min = -2;
fD_max = 2;



%% Initialize figure
h_fig = figure(11);
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
        
        
        %% Initialize plot
        subplot(max_f_case_number, max_D_case_number, count);
        hold on;
        [selected_bins_indices, selected_bins_centers] = get_selected_bins_indices();
        selected_bins_count = 1;
        % Use loaded values for the marginalized fD to set the x-limits
        x_lim_vec = [min(fD_mesh), max(fD_mesh)];
        
        
        %% Calculate the most probable values of f*D
        fD_values = zeros(1, lambda_count + 3);
        for lambda_ind = 1:lambda_count 
            bin = selected_bins_indices(lambda_ind);
            [mu_n, ~, ~, ~] = get_n_parameters(lambda_ind, bin);
            fD_values(lambda_ind) = mu_n * kBT / t_step;
        end;
        
        %% Calculate exact fD
        fD_values(lambda_count + 1: lambda_count + 3) = ...
            f_func(f_case_number, selected_bins_centers*L, L)...
            .* D_func(D_case_number, selected_bins_centers*L, L);
        
        %% Detect the longest interval and centrally increase its length
        factor = 100;
        fD_min = min(fD_values);
        fD_max = max(fD_values);
        fD_interval = fD_max - fD_min;
        fD_min = fD_min - fD_interval * (factor-1)/2;
        fD_max = fD_max + fD_interval * (factor-1)/2;
        
%         %% Manual adjustments
%         if f_case_number == 3 && D_case_number==2
%             val = 60;
%             fD_min = -val;
%             fD_max = val;
%         end;
        
        %% Prepare mesh
        fD_step = (fD_max - fD_min) / (fD_steps - 1);
        fD_mesh = fD_min:fD_step:fD_max;
        % Add special points
        fD_mesh = sort([fD_mesh, fD_values]);
        
        size1 = lambda_count;
        size2 = length(fD_mesh);
        fD_pdf_plot_data = zeros(size1, size2);
        % fD data
        for lambda_ind = 1:lambda_count
            bin = selected_bins_indices(lambda_ind);
            fD_pdf_plot_data(lambda_ind, :) = bin_fD_posterior_func (lambda_ind, bin, fD_mesh, t_step, kBT);
        end;


        %% Plot
        for lambda_ind = 1:lambda_count
            plot(fD_mesh, fD_pdf_plot_data(lambda_ind, :), '-',...
                'color', default_color_order(lambda_ind, :), 'LineWidth', 2);
        end;


        %% Simulated value
        y_lim_vec = [0, max(max(fD_pdf_plot_data)) * 1.1];
        plot(fD_values(lambda_count + 1) .* [1, 1], y_lim_vec, '--k', 'LineWidth', 2);
        ylim(y_lim_vec);    % I don't know why it changes scale without it
% 
% 
        %% Adjust
        if f_case_number == 1
            title(sprintf('D%i', D_case_number), 'FontWeight','Bold');
        end;
        if f_case_number == max_f_case_number
            xlabel('$\tilde f_j D_j$', 'interpreter', 'latex');
        end;
        if D_case_number == 1
            ylabel(sprintf('F%i', f_case_number), 'FontWeight','Bold');
        end;
%         ylabel('PDF');
%         x_lim_vec = xlim();
%         x_lim_vec(2) = x_lim_vec(2) * 1.1;
%         xlim(x_lim_vec);
%         x_lim_vec = [fD_min, fD_max];
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
output_filename = 'Disffusivity_times_force_posterior.pdf';
output_full_path = strcat(output_figures_folder, output_filename);
print(h_fig, output_full_path, '-dpdf', '-r0');








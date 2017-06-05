



%% Globals
global n_j;
global V_j;

%% Constants
height_to_width_ratio = 1;
marker_step = 10;
% % D_grad_steps = 1 + 2^10;
% % D_grad_min = -2;
% % D_grad_max = 2;


%% Initialize
[selected_bins_indices, selected_bins_centers] = get_selected_bins_indices();


%% Initialize figure
h_fig = figure(8);
clf;
set_article_figure_size(h_fig, height_to_width_ratio, 2);
default_color_order = get(gca,'colororder');
markers_list = {'-o','-s','-d','-^','-v'};


count = 1;
for f_case_number = 1:max_f_case_number
    for D_case_number = 1:max_D_case_number

        %% Get exact fD
        fD_val = f_func(f_case_number, selected_bins_centers*L, L) .* D_func(D_case_number, selected_bins_centers*L, L);
        
        %% Load data
        filename = sprintf('D_%i_f_%i_data.mat', D_case_number, f_case_number);
        full_path = strcat(output_data_folder, filename);
        load(full_path, '-mat');
        
        %% Initialize
        subplot(max_f_case_number, max_D_case_number, count);
        hold on;


        %% Plot
        for lambda_ind = 1:lambda_count
            plot(fD_mesh, fD_pdf_plot_data(lambda_ind, :), '-',...
                'color', default_color_order(lambda_ind, :), 'LineWidth', 2);
        end;


        %% Simulated value
        y_lim_vec = [0, max(max(fD_pdf_plot_data)) * 1.1];
        if y_lim_vec(2) == 0
            y_lim_vec(2) = 0.1;
        end;
        plot(fD_val(end) .* [1, 1], y_lim_vec, '--k', 'LineWidth', 2);
        ylim(y_lim_vec);    % I don't know why it changes scale without it


        %% Adjust
        if f_case_number == 1
            title(sprintf('D%i', D_case_number), 'FontWeight','Bold');
        end;
        if f_case_number == max_f_case_number
            xlabel('$f_jD_j$', 'interpreter', 'latex');
        end;
        if D_case_number == 1
            ylabel(sprintf('F%i', f_case_number), 'FontWeight','Bold');
        end;
%         ylabel('PDF');
%         x_lim_vec = xlim();
%         x_lim_vec(2) = x_lim_vec(2) * 1.1;
%         xlim(x_lim_vec);
        x_lim_vec = [fD_min, fD_max];
        xlim(x_lim_vec);
        box on;
        
        count = count + 1;
%         pause(0.01);
    end;
end;


%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize',[fig_pos(3), fig_pos(4)]);

% Set filename
output_filename = 'Disffusivity_times_force_marginalized_posterior.pdf';
output_full_path = strcat(output_figures_folder, output_filename);
print(h_fig, output_full_path, '-dpdf', '-r0');












%% Globals
% global n_j;
% global V_j;5

%% Constants
load_constants;
height_to_width_ratio = 1;
f_case_number = 1;
D_case_number = 1;
marker_step = 10;
D_steps = 1 + 2^10;
D_PRECISION = 1e-5;
D_max = 2;
interval_scaling_factor = 50;



%% Initialize figure
h_fig = figure(5);
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
        
        %% Prepare mesh
        % Calculate the most probable values of D
        D_values = zeros(1, lambda_count + 3);
        for lambda_ind = 1:lambda_count 
            bin = selected_bins_indices(lambda_ind);
            [mu_n, kappa_n, nu_n, sigma2_n] = get_n_parameters(lambda_ind, bin);
            D_values(lambda_ind) = nu_n * sigma2_n / (2 * t_step * (2 + nu_n));
        end;
        
        % Calculate exact D
        D_values(lambda_count + 1 : lambda_count + 3) = D_func(D_case_number, selected_bins_centers * L, L);

        %% Detect the longest interval and centrally increase its length
        D_min = min(D_values);
        D_max = max(D_values);
        D_interval = D_max - D_min;
        D_min = max([D_min - D_interval * (interval_scaling_factor-1)/2, D_PRECISION]);    % D>=0
        D_max = D_max + D_interval * (interval_scaling_factor-1)/2;
        
        %% Prepare mesh
        D_step = (D_max - D_min) / (D_steps - 1);
        D_mesh = D_min:D_step:D_max;
        % Add special points calculated earlier
        D_mesh = sort([D_mesh, D_values]);
        
        s1 = lambda_count;
        s2 = length(D_mesh);
        D_pdf_plot_data = zeros(s1, s2);
        % D data
        for lambda_ind = 1:lambda_count
            D_pdf_plot_data(lambda_ind, :) =  bin_D_posterior_func (lambda_ind,...
                selected_bins_indices(lambda_ind), D_mesh, t_step);
        end;
        
        
        %% Checking the PDF normalization
        pdf_norm = zeros(1, lambda_count);
        for lambda_ind = 1:lambda_count
            pdf_norm(lambda_ind) = trapz(D_mesh, D_pdf_plot_data(lambda_ind, :));
        end;
        fprintf('Norm for F %i/%i, D %i/%i:\n', f_case_number, max_f_case_number, D_case_number, max_D_case_number);
        disp(pdf_norm);


        %% Plot
        for lambda_ind = 1:lambda_count
            plot(D_mesh, D_pdf_plot_data(lambda_ind, :), '-',...
                'color', default_color_order(lambda_ind, :), 'LineWidth', 2);
        end;


        %% Simulated value
        y_lim_vec = [0, max(max(D_pdf_plot_data)) * 1.1];
        plot(D_values(end) .* [1, 1], y_lim_vec, '--k', 'LineWidth', 2);
        ylim(y_lim_vec);    % I don't know why it changes scale without it


        %% Adjust
        if f_case_number == 1
            title(sprintf('D%i', D_case_number), 'FontWeight','Bold');
        end;
        if f_case_number == max_f_case_number
            xlabel('D');
        end;
        if D_case_number == 1
            ylabel(sprintf('F%i', f_case_number), 'FontWeight','Bold');
        end;
%         ylabel('PDF');
%         x_lim_vec = xlim();
%         x_lim_vec(2) = x_lim_vec(2) * 1.1;
%         xlim(x_lim_vec);
        x_lim_vec = [D_min, D_max];
        xlim(x_lim_vec);
        box on;
        
        count = count + 1;
    end;
end;


%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize',[fig_pos(3), fig_pos(4)]);

% Set filename
output_filename = 'Disffusivity_posterior.pdf';
output_full_path = strcat(output_figures_folder, output_filename);
print(h_fig, output_full_path, '-dpdf', '-r0');











% This function performs inference of D in all bins for the selected f and
% D case


%% Globals


%% Constants
load_constants;
height_to_width_ratio = 1;
output_filename = 'Inferred_marginalized_force_times_diffusivity_profile.pdf';
fD_PRECISION = 1e-7;
fD_ABS_MAX = 100; % 8 * 2 * 10;
f_case_number = 1;
D_case_number = 1;
% l_ind = 1;


%% Initialize
h_fig = figure(14);
set_article_figure_size(h_fig, height_to_width_ratio, 2);
default_color_order = get(gca,'colororder');
clf;

%%%%%%%%%%%%

count = 1;
for f_case_number = 1:max_f_case_number
    for D_case_number = 1:max_D_case_number
        subplot(max_f_case_number, max_D_case_number, count);
        hold on;


        %% Load data for the given D and f cases
        filename = sprintf('D_%i_f_%i_data.mat', D_case_number, f_case_number);
        full_path = strcat(output_data_folder, filename);
        load(full_path, '-mat');


        %% Plot
        for l_ind = 1:lambda_count
            x_bins_number = x_bins_number_saved(l_ind);
            x_bins_centers = x_bins_centers_saved{l_ind};
            inferred_MAP_fD = zeros(1, x_bins_number);
            % Choose bin
            for bin = 2:x_bins_number
                fprintf('Processing D: %i/%i, f: %i/%i, lambda: %i/%i, bin: %i/%i\n', D_case_number,...
                    max_D_case_number, f_case_number, max_f_case_number,...
                    l_ind, lambda_count, bin, x_bins_number);
                % bin_fD_marginalized_posterior_func (l_ind, bin, fD)
                % Find maximum a posteriori for D in this bin
                function_to_minimze = @(fD) - bin_fD_marginalized_posterior_func (l_ind, bin, fD);
%                 optim_options = optimset('Display','iter', 'TolX', fD_PRECISION);
                optim_options = optimset('TolX', fD_PRECISION);
                fD_max = fminbnd(function_to_minimze, - fD_ABS_MAX, fD_ABS_MAX, optim_options);
                inferred_MAP_fD(bin) = fD_max;
            end;
            
            plot(x_bins_centers(2:end), inferred_MAP_fD(2:end), 'LineWidth', 2, 'color', default_color_order(l_ind, :))
            pause(0.01);
            %     plot(x_bins_centers_saved{l_ind}/L, elements_in_bins_count{l_ind}, markers_list{l_ind}, ...
        %     'LineWidth', 2, 'color', default_color_order(l_ind, :),...
        %     'MarkerIndices', 1:marker_step:length(x_bins_centers_saved{l_ind}),...
        %     'MarkerSize', marker_size, 'MarkerFaceColor', default_color_order(l_ind, :));
        end;
        
        
        %% Theoretical equilibrium distribution
        x_theor_mesh = x_min:0.005:x_max;
        fD_data = D_func(D_case_number, x_theor_mesh, L) .* f_func(f_case_number, x_theor_mesh, L);
        % Normalize to a PDF
        plot(x_theor_mesh/L, fD_data, '--k', 'LineWidth', 2);


        %% Label
        if f_case_number == 1
            title(sprintf('D%i', D_case_number), 'FontWeight','Bold');
        end;
        if f_case_number == max_f_case_number
            xlabel('x');
        end;
        if D_case_number == 1
            ylabel(sprintf('F%i', f_case_number), 'FontWeight','Bold');
        end;
        box on;
        xlim([-0.5, 0.5]);
        ylim([-fD_ABS_MAX, fD_ABS_MAX]);

        
        count = count + 1;
%         pause(0.01);
    end;
end;

hold off;


%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize',[fig_pos(3), fig_pos(4)]);

% Set filename
output_full_path = strcat(output_figures_folder, output_filename);
print(h_fig, output_full_path, '-dpdf', '-r0');

















%% Globals
% global elements_in_bins_count_saved;
global N;
load_constants;
% global x_bins_widths_saved;
% global x_bins_centers_saved;

%% Constants
height_to_width_ratio = 1;
marker_step = 10;
x_over_L_min = -1/2;
x_over_L_max = 1/2;
x_over_L_step = 0.01;


%% Initialize
if bl_short_trajectory
    local_min_points_in_bin = min_points_in_bin_short_trajectories;   % for normal calculations
    local_min_points_in_bin_step = -1;
    output_filename = 'Point_density.pdf';
else
    local_min_points_in_bin = 7e4;   % for long trajectories
    local_min_points_in_bin_step = -1000;
    output_filename = 'Point_density_long_trajectory.pdf';
end




%% Initialize figure
h_fig = figure(2);
clf;
set_article_figure_size(h_fig, height_to_width_ratio, 2);
default_color_order = get(gca,'colororder');
% markers_list = {'-o','-s','-d','-^','-v'};


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
            % Rebin data with larger bins to get smooth density
            % distribution
            for min_bin = local_min_points_in_bin:local_min_points_in_bin_step:1
                disp(min_bin);
                [~, x_bins_centers, x_bins_number, x_bins_widths,...
                elements_in_bins_count, ~] = select_bins_adaptive_mesh(x_lambda(l_ind, :),...
                dx_lambda(l_ind, :), min_bin);
                if x_bins_number >1     % Check if there are at least two points for the density plot
                    break;
                end;
            end;
            
            %x mesh
%             x_mesh = x_bins_centers;
            % PDF
            point_density = elements_in_bins_count;
%             bin_widths = x_bins_widths_saved{l_ind};
%             area = point_density * x_bins_widths;
            point_density = point_density ./ x_bins_widths' ./ N;
            
%             point_density = point_density / trapz(x_bins_centers_saved{l_ind}/L, point_density);
        %     trapz(x_bins_centers_saved{l_ind}/L, point_density)
            plot(x_bins_centers', point_density, ...
            'LineWidth', 2, 'color', default_color_order(l_ind, :));
            %     plot(x_bins_centers_saved{l_ind}/L, elements_in_bins_count{l_ind}, markers_list{l_ind}, ...
        %     'LineWidth', 2, 'color', default_color_order(l_ind, :),...
        %     'MarkerIndices', 1:marker_step:length(x_bins_centers_saved{l_ind}),...
        %     'MarkerSize', marker_size, 'MarkerFaceColor', default_color_order(l_ind, :));
        end;
        
        
        %% Theoretical equilibrium distribution
        [~, U_data] = f_func(f_case_number, x_bins_centers, L);
        rho_data = exp(- U_data / kBT);
        % Normalize to a PDF
        rho_data = L * rho_data / trapz(x_bins_centers, rho_data);
        plot(x_bins_centers/L, rho_data, '--k', 'LineWidth', 2);


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
output_full_path = strcat(output_figures_folder, output_filename);
print(h_fig, output_full_path, '-dpdf', '-r0');




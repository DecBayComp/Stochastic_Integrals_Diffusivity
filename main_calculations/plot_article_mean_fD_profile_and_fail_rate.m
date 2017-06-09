


function plot_article_mean_fD_profile_and_fail_rate(data_struct, fig_count, bl_save_figures)
%% === Plot mean fD profile for each inference type ===


%% Constants
load_constants;
sublabel_x = 0.03;
sublabel_y = 0.08;
x_lim_vec = [0, x_max];
y_lim_vec_FR = [-0.02, 1.02];
y_lim_vec_profile = [-0.11, 0.325];
output_filename = 'Force_profile_fail_rate.pdf';
% Subplot parameters
SH = 0.06;
SV = 0.1;
ML = 0.06;
MR = 0.005;
MT = 0.05;
MB = 0.08;
% Skip some bins
plot_every = 2;



%% Plot
% Initalize plot
h_fig = figure(fig_count);
set_article_figure_size(h_fig, 2, 2, 1);
clf;
% Initalize subplots
subaxis(2, lambda_types_count, 1, 'SH', SH, 'SV', SV, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);
for lambda_type = 1:lambda_types_count
    %% == Fail rate plot ==
    subaxis(2, lambda_types_count, lambda_types_count + lambda_type);
    hold on;
    str_legend = {};
    % Plot each inference conventions
    for convention = 1:conventions_count
        plot(data_struct.x_bins_centers(1:plot_every:end),  data_struct.UR_fD(lambda_type, 1:plot_every:end, convention),...
            strcat('-', markers_list{convention}), 'LineWidth', line_width, 'color', color_sequence(convention, :), 'markers', marker_size);
        % Update legend
        str_legend{end + 1} = conventions_names{convention};
    end
    % Confidence level (theory)
    plot(x_lim_vec, x_lim_vec * 0 + (1 - CONF_LEVEL), 'k--', 'LineWidth', line_width);

    %% Adjust
    xlim(x_lim_vec);
    ylim(y_lim_vec_FR);
    box on;
    
    xlabel('$x$', 'interpreter', 'latex');
    ylabel('Fail rate', 'interpreter', 'latex');
    % Subplot label
    text(sublabel_x, sublabel_y, strcat('(', char('e' + lambda_type - 1), ')'), 'Units', 'Normalized', 'VerticalAlignment', 'Top');
    
    % Legend
    if lambda_type == 2
        legend(str_legend, 'location', 'northwest', 'FontSize', legend_font_size);
        legend boxon;
    end;   
 
    
    %% == Profile plot ==
    subaxis(2, lambda_types_count, lambda_type);
    hold on;
    % Plot each convention
    for convention = 1:conventions_count
        plot(data_struct.x_bins_centers(1:plot_every:end),  data_struct.MAP_fD_mean(lambda_type, 1:plot_every:end, convention, 1),...
            strcat('-', markers_list{convention}), 'color', color_sequence(convention, :), 'LineWidth', line_width, 'markers', marker_size);
    end;
    % True profile (theory)
    h_theor = plot(data_struct.x_fine_mesh, data_struct.fD_theor_fine_data, '--k', 'LineWidth', line_width);

    %% Adjust
    xlim(x_lim_vec);
    ylim(y_lim_vec_profile);
    box on;
    
    xlabel('$x$', 'interpreter', 'latex');
    ylabel('$fD$', 'interpreter', 'latex');
    str_title = {'$\lambda^* = 0$', '$\lambda^* = 0.5$', '$\lambda^* = 1$', 'Random $\lambda^*$'};
    title(str_title{lambda_type}, 'interpreter', 'latex');
    % Subplot label
    text(sublabel_x, sublabel_y, strcat('(', char('a' + lambda_type - 1), ')'), 'Units', 'Normalized', 'VerticalAlignment', 'Top');
    
    % Push theoretical curve back
    uistack(h_theor, 'bottom');
end;


%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize', [fig_pos(3), fig_pos(4)]);
% Set filename
output_full_path = strcat(output_figures_folder, output_filename);
if bl_save_figures
    print(h_fig, output_full_path, '-dpdf', '-r0');
end;










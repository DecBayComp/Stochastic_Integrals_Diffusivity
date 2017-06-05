


function fig_count = plot_article_mean_fD_profile(data_struct, fig_count, bl_save_figures)


%% Constants
load_constants;
sublabel_x = 0.02;
sublabel_y = 0.99;


%% === Plot mean fD profile for each inference type ===
fig_count = fig_count + 1; 
h_fig = figure(fig_count);
set_article_figure_size(h_fig, 2, 2, 1);
clf;
x_lim_vec = [0, x_max];
SH = 0.08;
SV = 0.1;
ML = 0.07;
MR = 0.005;
MT = 0.05;
MB = 0.08;
for l_ind = 1:lambda_count
    %% Success rate plot
    h_sub = subaxis(2, lambda_count, lambda_count + l_ind, 'SH', SH, 'SV', SV, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);
    hold on;
    % Divine
    plot(data_struct.x_bins_centers_saved{l_ind},  data_struct.UR_fwd_divine{l_ind},...
        strcat('-', markers_list{1}), 'LineWidth', line_width);
    % Ito
    plot(data_struct.x_bins_centers_saved{l_ind},  data_struct.UR_Ito{l_ind},...
        strcat('-', markers_list{2}), 'LineWidth', line_width);
    % Stratonovich
    plot(data_struct.x_bins_centers_saved{l_ind},  data_struct.UR_fwd_Stratonovich{l_ind},...
        strcat('-', markers_list{3}), 'LineWidth', line_width);
    % Marginalized
    plot(data_struct.x_bins_centers_saved{l_ind},  data_struct.UR_fwd_marginalized{l_ind},...
        'LineWidth', line_width);
%     % Direct Hanggi
%     plot(data_struct.x_bins_centers_saved{l_ind},  data_struct.UR_bck_Hanggi{l_ind},...
%         'LineWidth', line_width);
    % Adjust
    xlim(x_lim_vec);
    ylim([-0.02, 1.02]);
    xlabel('x');
    ylabel('Fail rate');
    box on;
    pos = get(gca, 'Position');
    % Add subplot label
    text(sublabel_x, sublabel_y, strcat('(', char('d' + l_ind - 1), ')'), 'Units', 'Normalized', 'VerticalAlignment', 'Top');
%     set(gca, 'LooseInset', get(gca,'TightInset'))
    if l_ind == 2
%         legend_position = [0.411, 0.583, 0.067, 0.127];
        str_legend = {'Orcl', 'Ito', 'Str', 'Mar'};
        h_leg = legend(str_legend, 'location', 'southwest');
        legend boxon;
    end;   
 
    
    %% Profile plot
    h_sub = subaxis(2, lambda_count, l_ind);
%     subplot(2, lambda_count, lambda_count + l_ind);
    hold on;
    % Divine
    plot(data_struct.x_bins_centers_saved{l_ind},  data_struct.MAP_fwd_fD_divine_mean{l_ind}(1, :),...
        strcat('-', markers_list{1}), 'LineWidth', line_width);
    % Ito
    plot(data_struct.x_bins_centers_saved{l_ind},  data_struct.MAP_fD_Ito_mean{l_ind}(1, :),...
        strcat('-', markers_list{2}), 'LineWidth', line_width);
    % Stratonovich
    plot(data_struct.x_bins_centers_saved{l_ind},  data_struct.MAP_fwd_fD_Stratonovich_mean{l_ind}(1, :),...
        strcat('-', markers_list{3}), 'LineWidth', line_width);
    % Marginalized
    plot(data_struct.x_bins_centers_saved{l_ind},  data_struct.MAP_fwd_fD_marginalized_mean{l_ind}(1, :),...
        'LineWidth', line_width);
    % Theory
    plot(data_struct.x_fine_mesh, data_struct.fD_theor_fine_data, '--k', 'LineWidth', line_width);
%     % Direct Hanggi
%     plot(data_struct.x_bins_centers_saved{l_ind},  -data_struct.MAP_bck_fD_Hanggi_mean{l_ind}(1, :),...
%         'LineWidth', line_width);
    % Adjust
    xlim(x_lim_vec);
%     ylim([-0.1, 0.3]);
    ylim([-0.07, 0.275]);
    xlabel('x');
    ylabel('fD');
    str_title = {'$\lambda^* = 0$', '$\lambda^* = 1/2$', '$\lambda^* = 1$'};
    title(str_title{l_ind}, 'interpreter', 'latex');
    box on;
    % Add subplot label
    text(sublabel_x, sublabel_y, strcat('(', char('a' + l_ind - 1), ')'), 'Units', 'Normalized', 'VerticalAlignment', 'Top');
%     set(gca, 'LooseInset', get(gca,'TightInset'))
        % Legend
    

end;

%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize', [fig_pos(3), fig_pos(4)]);
% Set filename
output_filename = 'Force_profile_fail_rate.pdf';
output_full_path = strcat(output_figures_folder, output_filename);
if bl_save_figures
    print(h_fig, output_full_path, '-dpdf', '-r0');
end;
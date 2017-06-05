



function fig_count = plot_article_fD_error_cones(data_struct, fig_count, bl_save_figures)


%% Constants
load_constants;


%% === Plot fD error cones ===
fig_count = fig_count + 1; 
h_fig = figure(fig_count);
set_article_figure_size(h_fig, 1, 2, 1);
clf;
x_lim_vec = [-1, 1] * 0.165;
y_lim_vec = [0, 1.5];
spacing = 0.08;
ML = 0.07;
MR = 0.02;
MT = 0.08;
MB = 0.17;
for l_ind = 1:lambda_count
    h_sub = subaxis(1, lambda_count, l_ind, 'Spacing', spacing, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);
    hold on;
    % Filter indices from one best period
    w = 10.0;
    x_left = (1/2 + 2*1)/w;
    x_right = (1/2 + 2*2)/w;
    indices = data_struct.x_bins_centers_saved{l_ind} >= x_left & data_struct.x_bins_centers_saved{l_ind} <= x_right;
%     % Divine
%     difference = (data_struct.MAP_fwd_fD_divine_mean{l_ind}(1, indices) - data_struct.fD_theor_data{l_ind}(indices)') / kBT;
%     half_error = (data_struct.MAP_fwd_fD_divine_mean{l_ind}(2, indices) + data_struct.MAP_fwd_fD_divine_mean{l_ind}(3, indices))/2 / kBT;
%     plot(half_error, difference, markers_list{1}, 'markers', marker_size);
%     % Ito
%     difference = (data_struct.MAP_fD_Ito_mean{l_ind}(1, indices) - data_struct.fD_theor_data{l_ind}(indices)') / kBT;
%     half_error = (data_struct.MAP_fD_Ito_mean{l_ind}(2, indices) + data_struct.MAP_fD_Ito_mean{l_ind}(3, indices))/2 / kBT; 
%     plot(difference, half_error, markers_list{2}, 'markers', marker_size);
    % Stratonovich
    difference = (data_struct.MAP_fwd_fD_Stratonovich_mean{l_ind}(1, indices) - data_struct.fD_theor_data{l_ind}(indices)') / kBT;
    half_error = (data_struct.MAP_fwd_fD_Stratonovich_mean{l_ind}(2, indices) + data_struct.MAP_fwd_fD_Stratonovich_mean{l_ind}(3, indices))/2 / kBT;
    difference_rel = abs(difference ./ half_error);
    plot(data_struct.x_bins_centers_saved{l_ind}(indices), difference_rel, markers_list{3}, 'markers', marker_size);
    1;
%     % Marginalized
%     difference = (data_struct.MAP_fwd_fD_marginalized_mean{l_ind}(1, indices) - data_struct.fD_theor_data{l_ind}(indices)') / kBT;
%     half_error = (data_struct.MAP_fwd_fD_marginalized_mean{l_ind}(2, indices) + data_struct.MAP_fwd_fD_marginalized_mean{l_ind}(3, indices))/2 / kBT;
%     plot(difference, half_error, markers_list{4}, 'markers', marker_size);
    
%     % Identity
%     x_temp = [ - y_lim_vec(2), 0, y_lim_vec(2)];
%     y_temp = [-1, 0, 1] .* x_temp;
%     plot(x_temp, y_temp, '--k');
    
%      % Legend
%     if l_ind == 1
%         str_legend = {'Orcl', 'Ito', 'Str', 'Mar'};
%         legend(str_legend, 'location', 'southwest');
%     end;
    
    % Adjust
%     xlim(x_lim_vec);
    ylim(y_lim_vec);
%     ylabel('fD bias');
%     xlabel('Half confidence interval');
    title(sprintf('$\\lambda^* = %.2f$', lambda_array(l_ind)), 'interpreter', 'latex');
end;


%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize', [fig_pos(3), fig_pos(4)]);
% Set filename
output_filename = 'Force_error_cones.pdf';
output_full_path = strcat(output_figures_folder, output_filename);
if bl_save_figures
    print(h_fig, output_full_path, '-dpdf', '-r0');
end;
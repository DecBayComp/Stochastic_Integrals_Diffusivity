

%% Constants
color_sequence = {'b', 'r', 'k'};


%% Calculating theoretical D and f
x_plot_points_number = 100;
x_plot_step = (x_max-x_min)/(x_plot_points_number -1);
x_plot_mesh = x_min:x_plot_step:x_max;
% alpha_map = alpha_func(x_plot_mesh);
D_map_theory = D_func(x_plot_mesh);
f_map_theory = f_func(x_plot_mesh);


%% Plotting D

figure(3);
clf;
hold all;

legend_string = cell(1, lambda_count + 1);

% D
% Plotting theory
plot(x_plot_mesh/L, D_map_theory, 'LineWidth', 2);
legend_string{1} = 'theory';
% Plotting inference
for l_ind = 1:lambda_count
    plot(x_bins_centers_saved{l_ind}/L, D_map_inferred{l_ind}, 'LineWidth', 2);
    legend_string{l_ind + 1} = sprintf('\\lambda = %.2f', lambda_array(l_ind));
end;

xlabel('x/L');
ylabel('D');
ylim_vec = ylim();
ylim_vec(1) = 0;
ylim(ylim_vec);
legend(legend_string, 'Location', 'best');
hold off;



%% Plotting f
% Saving legend position
if exist('h_force_legend', 'var')
%     figure(fig_force);
%     h_force_legend = legend;
    legend_position = h_force_legend.Position;
end;
fig_force = figure(4);
clf;
hold all;
legend_string = cell(1, 2 * lambda_count + 1);

% Plotting theory
plot(x_plot_mesh/L, f_map_theory, 'm', 'LineWidth', 2);
legend_string{1} = 'theory';
legends_saved = 1;
% Plotting inference
for l_ind = 1:lambda_count
    % Plotting the force with lambda = 0
    plot(x_bins_centers_saved{l_ind}/L, f_Pb2_map_inferred{l_ind}, color_sequence{l_ind}, 'LineWidth', 2);
    legends_saved = legends_saved + 1;
    legend_string{legends_saved} = sprintf('\\lambda = %.2f, max', lambda_array(l_ind));
    % Plotting the force with lambda = 1
    plot(x_bins_centers_saved{l_ind}/L, f_Pb2_map_inferred{l_ind} - f_lambda_map_inferred{l_ind},...
        strcat(color_sequence{l_ind}, '--'), 'LineWidth', 2);
    legends_saved = legends_saved + 1;
    legend_string{legends_saved} = sprintf('\\lambda = %.2f, min', lambda_array(l_ind));
end;

xlabel('x/L');
ylabel('a(x)');
ylim_vec = ylim();
% % ylim_vec(1) = 0;
% % ylim(ylim_vec);
if exist('legend_position', 'var')
    h_force_legend = legend(legend_string, 'Position', legend_position);
else
    h_force_legend = legend(legend_string, 'Location', 'best');
end;
hold off;

















%% Calculating D
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
% % % % Alpha
% % % subplot(1, 2, 1);
% % % plot(x_plot_mesh/L, alpha_map, 'LineWidth', 2);
% % % xlabel('x/L');
% % % ylabel('\alpha');
% % % ylim_vec = ylim();
% % % ylim_vec(1) = 0;
% % % ylim(ylim_vec);

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

figure(4);
clf;
hold all;
% % % % Alpha
% % % subplot(1, 2, 1);
% % % plot(x_plot_mesh/L, alpha_map, 'LineWidth', 2);
% % % xlabel('x/L');
% % % ylabel('\alpha');
% % % ylim_vec = ylim();
% % % ylim_vec(1) = 0;
% % % ylim(ylim_vec);

% % legend_string = cell(1, lambda_count + 1);

% f
% Plotting theory
plot(x_plot_mesh/L, f_map_theory, 'LineWidth', 2);
% % legend_string{1} = 'theory';
% Plotting inference
for l_ind = 1:lambda_count
    plot(x_bins_centers_saved{l_ind}/L, f_map_inferred{l_ind}, 'LineWidth', 2);
% %     legend_string{l_ind + 1} = sprintf('\\lambda = %.2f', lambda_array(l_ind));
end;

xlabel('x/L');
ylabel('a(x)');
ylim_vec = ylim();
% % ylim_vec(1) = 0;
% % ylim(ylim_vec);
legend(legend_string, 'Location', 'best');
hold off;













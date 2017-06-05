


%% Constants
x_over_L_min = -1/2;
x_over_L_max = 1/2;
x_over_L_step = 0.1;
height_to_width_ratio = 1/1.5;
f_case_number = 1;
D_case_number = 1;
lambda_ind = 1;
plot_step = 1;




%% Initialize figure
h_fig = figure(1);
set_article_figure_size(h_fig, height_to_width_ratio, 3);


%% Load data for the given D and f cases
filename = sprintf('D_%i_f_%i_data.mat', D_case_number, f_case_number);
full_path = strcat(output_data_folder, filename);
load(full_path, '-mat');


%% Plot trajectory
plot(t_mesh(1:plot_step:end-1), x_lambda(lambda_ind, 1:plot_step:end), 'LineWidth', 2);
xlabel('t');
ylabel('x(t)');
x_lim_vec = [0, t_mesh(end)];
xlim(x_lim_vec);
% title(sprintf('Simulated trajectory for F%i, D%i, \\lambda = %.1f', f_case_number,...
%     D_case_number, (lambda_ind-1) * 0.5));


%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize',[fig_pos(3), fig_pos(4)]);

% Set filename
output_filename = 'Sample_trajectory.pdf';
output_full_path = strcat(output_figures_folder, output_filename);
print(h_fig, output_full_path, '-dpdf', '-r0');






















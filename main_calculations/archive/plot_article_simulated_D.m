


%% Globals
global L;
global max_D_case_number;

%% Constants
marker_step = 20;
x_over_L_min = -1/2;
x_over_L_max = 1/2;
x_over_L_step = 0.1/marker_step;
selected_x_over_L = 1/4;
golden_ratio = (1+sqrt(5))/2;
height_to_width_ratio = 1/1.5;


%% Initialize figure
h_fig = figure(10);
set_article_figure_size(h_fig, height_to_width_ratio, 3);
% % Set legend position
% if exist('h_legend', 'var')
%     legend_position = h_legend.Position;
% end;
legend_position = [0.3297    0.6677    0.1448    0.1820];


clf;
hold on;
default_color_order = get(gca,'colororder');
markers_list = {'-o','-s','-d','-^','-v'};


%% Initialize x mesh
x_mesh = (x_over_L_min:x_over_L_step:x_over_L_max) * L;
x_mesh_length = length(x_mesh);

%% Calculate force for each case
D_data = zeros(max_D_case_number, x_mesh_length);
for D_case = 1:max_D_case_number
    D_data(D_case, :) = D_func(D_case, x_mesh, L);
end;


%% Plot
for D_case = 1:max_D_case_number
    plot(x_mesh/L, D_data(D_case, :), markers_list{D_case}, 'LineWidth', 2,...
    'color', default_color_order(D_case, :), ...
    'MarkerSize', marker_size, 'MarkerFaceColor', default_color_order(D_case, :),...
    'MarkerIndices', 1:marker_step:x_mesh_length);
end;


%% Label
y_lim_vec = [1, 2];
ylim(y_lim_vec);
box on;
grid on;
xlabel('x');
ylabel('D(x)');


%% Add a vertical line at a selected x location
x = selected_x_over_L * ones(1, 2);
y = y_lim_vec;
plot(x, y, '--k', 'LineWidth', 2);


%% Legend
str_legend = {'D1', 'D2', 'D3', 'D4'};
h_legend = legend(str_legend, 'FontSize', font_size - 5, 'location', 'northeastoutside'); % 'Position', legend_position);
% h_legend = legend(str_legend, 'location', 'northwest');
% % Keep legend in place
% if exist('h_legend', 'var')
%     set(h_legend, 'Position', legend_position);
% end


%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize',[fig_pos(3), fig_pos(4)]);

% Set filename
output_filename = 'Simulated_diffusivity.pdf';
output_full_path = strcat(output_figures_folder, output_filename);
print(h_fig, output_full_path, '-dpdf', '-r0');








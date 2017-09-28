%% Plot the simulated diffusivity and force profiles



%% Constants
load_constants;
load_color_scheme;
x_marker_step = 0.05 * L;
marker_step = 10;
x_step = x_marker_step/marker_step;
output_D_filename = 'Simulated_D_a_b.pdf';

% Subplot parameters
spacing = 0.08;
ML = 0.07;
MR = 0.0125;
MT = 0.08;
MB = 0.17;



%% Create mesh
x_mesh = x_min:x_step:x_max;
x_mesh_length = length(x_mesh);

%% Calculate
% D
D_data = D_func(selected_D_case, x_mesh, L);
% Force f
f_data = f_func(selected_f_case, x_mesh, L);
% Local drag
a_data = f_data / gamma_drag;	% in um/s
% Diffusivity b(x)
b_data = sqrt(2*D_data);



%% Initialize figure
% Common figure
h_fig = figure(10);
clf;
set_article_figure_size(h_a_fig, 1, 2, 1);

% Initialize subplots
h_sub = subaxis(1, 3, 1, 'Spacing', spacing, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);



%% a figure
% Initialize subplot
h_sub = subaxis(1, 3, 1);

% Plot
% color = 
plot(x_mesh, a_data, '-o', 'LineWidth', 2,...
    'color', standard_colors(1).DeepBlue, ...
	'MarkerSize', marker_size, 'MarkerFaceColor', standard_colors(1).DeepBlue,...
    'MarkerIndices', 1:marker_step:x_mesh_length);

% Label
ylim([0, max(a_data)]);
box on;
grid on;
xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
ylabel('$a$, $\mu \mathrm{m/s}$', 'interpreter', 'latex');



%% Diffusivity figure
% Initialize sbuplot
h_sub = subaxis(1, 3, 2);

% Plot
% color = 
plot(x_mesh, D_data, '-o', 'LineWidth', 2,...
    'color', standard_colors(1).DeepBlue, ...
	'MarkerSize', marker_size, 'MarkerFaceColor', standard_colors(1).DeepBlue,...
    'MarkerIndices', 1:marker_step:x_mesh_length);

% Label
box on;
grid on;
xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
ylabel('$D$, $\mu \mathrm{m^2/s}$', 'interpreter', 'latex');



%% b figure
% Initialize sbuplot
h_sub = subaxis(1, 3, 3);

% Plot
% color = 
plot(x_mesh, b_data, '-o', 'LineWidth', 2,...
    'color', standard_colors(1).DeepBlue, ...
	'MarkerSize', marker_size, 'MarkerFaceColor', standard_colors(1).DeepBlue,...
    'MarkerIndices', 1:marker_step:x_mesh_length);

% Label
y_lim_vec = [min(b_data), max(b_data)];
ylim(y_lim_vec);
box on;
grid on;
xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
ylabel('$b$, $\mu \mathrm{m \cdot s^{-1/2}}$', 'interpreter', 'latex');




%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize', [fig_pos(3), fig_pos(4)]);
% Set filename
output_full_path = strcat(output_figures_folder, output_D_filename);
% Print
print(h_fig, output_full_path, '-dpdf', '-r0');







%% Plot the simulated diffusivity and force profiles



function plot_simulated_diffusivity_and_force(fig_count, bl_save_figures)



%% Constants
load_constants;
load_color_scheme;
x_marker_step = 0.05 * L;
marker_step = 10;
x_step = x_marker_step/marker_step;
output_D_filename = 'simulated_D_a_b.pdf';

% Subplot parameters
rows = 1;
cols = 4;
spacing = 0.08;
ML = 0.07;
MR = 0.0125;
MT = 0.09;
MB = 0.17;

% Label params
sublabel_x = -0.03;
sublabel_y = 1.12;

x_tick_increment = 0.25;



%% Create mesh
x_mesh = x_min:x_step:x_max;
x_mesh_length = length(x_mesh);

%% Calculate
% D
D_data = D_func(selected_D_case, x_mesh, L);
% Force f
f_data_force = f_func(enum_force_case, x_mesh, L);
f_data_no_force = f_func(enum_no_force_case, x_mesh, L);
% Local drag
a_data_force = f_data_force / gamma_drag;	% in um/s
a_data_no_force = f_data_no_force / gamma_drag;	% in um/s
% Diffusivity b(x)
b_data = sqrt(2*D_data);

a_y_lim_vec = max(abs([a_data_force, a_data_no_force])) * [-1, 1];



%% Initialize figure
% Common figure
h_fig = figure(fig_count);
clf;
set_article_figure_size(h_fig, 1, 2, 1);

% Initialize subplots
h_sub = subaxis(rows, cols, 1, 'Spacing', spacing, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);



%% 'a' figure for FORCE case
% Initialize subplot
h_sub = subaxis(1);

% Plot
% color = 
plot(x_mesh, a_data_force, '-o', 'LineWidth', 2,...
    'color', standard_colors(1).DeepBlue, ...
	'MarkerSize', marker_size, 'MarkerFaceColor', standard_colors(1).DeepBlue,...
    'MarkerIndices', 1:marker_step:x_mesh_length);

% Adjust
% max_a = max(a_data_force);
% if max_a ~= 0
% 	ylim([0, max_a]);
% else
% 	ylim('auto');
% end;
ylim(a_y_lim_vec);
box on;
grid on;
xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
ylabel('$a$, $\mu \mathrm{m/s}$', 'interpreter', 'latex');

title('Local-force mod.', 'interpreter', 'latex');

% Modify ticks
set(gca,'xtick', x_min:x_tick_increment:x_max);

% Subplot label
text(sublabel_x, sublabel_y, 'A', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', subplot_label_font_size);



%% 'a' figure for NO FORCE case
% Initialize subplot
h_sub = subaxis(2);

% Plot
% color = 
plot(x_mesh, a_data_no_force, '-o', 'LineWidth', 2,...
    'color', standard_colors(1).DeepBlue, ...
	'MarkerSize', marker_size, 'MarkerFaceColor', standard_colors(1).DeepBlue,...
    'MarkerIndices', 1:marker_step:x_mesh_length);

% Adjust
% max_a = max(a_data_force);
% if max_a ~= 0
% 	ylim([0, max_a]);
% else
% 	ylim('auto');
% end;
ylim(a_y_lim_vec);
box on;
grid on;
xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
ylabel('$a$, $\mu \mathrm{m/s}$', 'interpreter', 'latex');

title('Spurious-force mod.', 'interpreter', 'latex');

% Modify ticks
set(gca,'xtick', x_min:x_tick_increment:x_max);

% Subplot label
text(sublabel_x, sublabel_y, 'B', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', subplot_label_font_size);



%% Diffusivity figure
% Initialize sbuplot
h_sub = subaxis(3);

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

% Modify ticks
set(gca,'xtick', x_min:x_tick_increment:x_max);

% Subplot label
text(sublabel_x, sublabel_y, 'C', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', subplot_label_font_size);



%% b figure
% Initialize sbuplot
h_sub = subaxis(4);

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

% Modify ticks
set(gca,'xtick', x_min:x_tick_increment:x_max);

% Subplot label
text(sublabel_x, sublabel_y, 'D', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', subplot_label_font_size);




%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize', [fig_pos(3), fig_pos(4)]);
% Set filename
output_full_path = strcat(output_figures_folder, output_D_filename);
% Print
if bl_save_figures
	print(h_fig, output_full_path, '-dpdf', '-r0');
end;







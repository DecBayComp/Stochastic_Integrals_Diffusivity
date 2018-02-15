%% Plot the simulated diffusivity and force profiles



function plot_simulated_diffusivity_and_force(fig_count, bl_save_figures)



%% Constants
load_constants;
load_color_scheme;
x_step = L * 1e-3;
x_marker_step = L;
marker_step = round(x_marker_step/x_step);
output_D_filename = 'simulated_D_a_b.pdf';

% Subplot parameters
rows = 2;
cols = 2;
spacing = 0.08;
ML = 0.19;
MR = 0.03;
MT = 0.08;
MB = 0.16;
SH = 0.2;
SV = 0.19;

% Label params
sublabel_x = 0;
sublabel_y = 1.3;

% x_tick_increment = 0.25;



%% Create mesh
x_mesh = x_min:x_step:x_max;
x_mesh_length = length(x_mesh);

%% Calculate
% D
[D_data, D_der_data] = D_func(selected_D_case, x_mesh, L);

% Force f
f_data_force = f_func(enum_force_case, x_mesh, L);
f_data_no_force = f_func(enum_no_force_case, x_mesh, L);

% Conservative force [3 x mesh_size]
a_data_force = [1; 1; 1] * f_data_force / gamma_drag;	% in um/s
a_data_no_force = [1; 1; 1] * f_data_no_force / gamma_drag;	% in um/s

% Spurious force [3 x mesh_size]
spur_data = [0; 0.5; 1] * D_der_data;

% Total force
alpha_data_no_force = spur_data + a_data_no_force;
alpha_data_force = spur_data + a_data_force;


% Diffusivity b(x)
b_data = sqrt(2*D_data);

a_y_lim_vec = max(max(abs([alpha_data_force, alpha_data_no_force]))) * [-1, 1];



%% Initialize figure
% Common figure
h_fig = figure(fig_count);
clf;
set_article_figure_size(h_fig, 2, 1, 0.6);

% Initialize subplots
h_sub = subaxis(rows, cols, 1, 'SH', SH, 'SV', SV, 'ML', ML, 'MR', MR, 'MT', MT, 'MB', MB);



%% Diffusivity figure
% Initialize sbuplot
h_sub = subaxis(1);

% Plot
% color = 
plot(x_mesh, D_data, '-', 'LineWidth', 2,...
    'color', standard_colors(1).DeepBlue); %, ...
	%'MarkerSize', marker_size, 'MarkerFaceColor', standard_colors(1).DeepBlue,...
    % 'MarkerIndices', 1:marker_step:x_mesh_length);

% Label
box on;
% grid on;
% xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
ylabel('$D$, $\mu \mathrm{m^2/s}$', 'interpreter', 'latex');

% Modify ticks
set(gca, 'FontSize', font_size);

% set(gca,'xtick', x_min:x_tick_increment:x_max);

% Subplot label
text(sublabel_x, sublabel_y, 'A', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', subplot_label_font_size);



%% b figure
% Initialize sbuplot
h_sub = subaxis(2);

% Plot
% color = 
plot(x_mesh, b_data, '-', 'LineWidth', 2,...
    'color', standard_colors(1).DeepBlue, ...
	'MarkerSize', marker_size, 'MarkerFaceColor', standard_colors(1).DeepBlue,...
    'MarkerIndices', 1:marker_step:x_mesh_length);

% Label
y_lim_vec = [min(b_data), max(b_data)];
ylim(y_lim_vec);
box on;
% grid on;
% xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
ylabel('$b$, $\mu \mathrm{m \cdot s^{-1/2}}$', 'interpreter', 'latex');

% Modify ticks
set(gca, 'FontSize', font_size);
% set(gca,'xtick', x_min:x_tick_increment:x_max);

% Subplot label
text(sublabel_x, sublabel_y, 'B', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', subplot_label_font_size);



%% 'alpha' for NO FORCE and FORCE
force_data_array = {alpha_data_no_force, alpha_data_force};
title_array = {'No force', 'Force'};
% Initialize subplot
for f_ind = 1:2
h_sub = subaxis(2 + f_ind);

hold on;

% Parse simulated fixed-lambda^*
for convention = 1:3
    plot(x_mesh, force_data_array{f_ind}(convention, :), strcat('-', markers_list{convention}),...
        'color', color_sequence(convention, :), 'MarkerIndices', 1:marker_step:x_mesh_length, 'markers', marker_size-2, 'linewidth', line_width);
end

% Adjust
ylim(a_y_lim_vec);
box on;
% grid on;
axis manual;
xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
ylabel('$\alpha$, $\mu \mathrm{m/s}$', 'interpreter', 'latex');

title(title_array{f_ind}, 'interpreter', 'latex');

% Modify ticks
set(gca, 'FontSize', font_size);
% set(gca,'xtick', x_min:x_tick_increment:x_max);

% Zero level
h_zero = plot(xlim(), [0,0], '-', 'linewidth', line_width_theor, 'color', axes_color);
uistack(h_zero, 'bottom');

% Subplot label
text(sublabel_x, sublabel_y, char('B' + f_ind), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', subplot_label_font_size);

end



% % % % % 
% % % % % 
% % % % % 
% % % % % %% 'a' figure for FORCE case
% % % % % % Initialize subplot
% % % % % h_sub = subaxis(4);
% % % % % 
% % % % % % Plot
% % % % % % color = 
% % % % % plot(x_mesh, a_data_force, '-o', 'LineWidth', 2,...
% % % % %     'color', standard_colors(1).DeepBlue, ...
% % % % % 	'MarkerSize', marker_size, 'MarkerFaceColor', standard_colors(1).DeepBlue,...
% % % % %     'MarkerIndices', 1:marker_step:x_mesh_length);
% % % % % 
% % % % % % Adjust
% % % % % % max_a = max(a_data_force);
% % % % % % if max_a ~= 0
% % % % % % 	ylim([0, max_a]);
% % % % % % else
% % % % % % 	ylim('auto');
% % % % % % end;
% % % % % ylim(a_y_lim_vec);
% % % % % box on;
% % % % % grid on;
% % % % % xlabel('$x$, $\mu \mathrm{m}$', 'interpreter', 'latex');
% % % % % ylabel('$a$, $\mu \mathrm{m/s}$', 'interpreter', 'latex');
% % % % % 
% % % % % % title('Local-force mod.', 'interpreter', 'latex');
% % % % % 
% % % % % % Modify ticks
% % % % % set(gca, 'FontSize', font_size);
% % % % % % set(gca,'xtick', x_min:x_tick_increment:x_max);
% % % % % 
% % % % % % Subplot label
% % % % % text(sublabel_x, sublabel_y, 'D', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', subplot_label_font_size);




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
end







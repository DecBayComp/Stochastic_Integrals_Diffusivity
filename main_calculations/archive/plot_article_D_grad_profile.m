


function fig_count = plot_article_D_grad_profile(data_struct, fig_count, bl_save_figures)


%% Constants
load_constants;

%% === Plot D' average regularized profile ===
%% Calculate simple D' with MAP Dvalues for lambda = 1
% lambda_type = enum_lambda_Stratonovich;
x_lim_vec = [0, x_max];
y_lim_vec = [-1.1, 0.65];
x_bins_steps = data_struct.x_bins_centers(2:end) - data_struct.x_bins_centers(1:end - 1);
simple_D_grad = (data_struct.MAP_D(2:end, 1) - data_struct.MAP_D(1:end-1, 1)) ./ x_bins_steps;
x_grad_mesh = data_struct.x_grad_mesh;

fig_count = fig_count + 1; 
h_fig = figure(fig_count);
set_article_figure_size(h_fig, 1, 1, 1);
clf;
hold on;
% Simple difference
plot(x_grad_mesh, simple_D_grad, markers_list{1}, 'LineWidth', line_width, 'markers', marker_size-2);
% Regularized gradient
plot(x_grad_mesh, data_struct.MAP_D_grad_regular,  markers_list{2}, 'LineWidth', line_width, 'markers', marker_size-2);
% Regularized interpolated gradient
plot(data_struct.x_bins_centers, data_struct.MAP_D_grad_regular_interp, 'LineWidth', line_width, 'markers', marker_size-2);
% Theory
plot(data_struct.x_fine_mesh, data_struct.D_grad_theor_fine_data, 'k--', 'LineWidth', 2);
% Adjust
xlabel('x', 'interpreter', 'latex');
ylabel('$\nabla D$', 'interpreter', 'latex');
xlim(x_lim_vec);
ylim(y_lim_vec);
box on;
title(sprintf('D gradient profile for $\\lambda^* = %.2f$', data_struct.lambda), 'interpreter', 'latex');
% Legend
str_legend_local = {'FDif', 'Reg', 'RegI'};
legend(str_legend_local, 'location', 'south', 'interpreter', 'latex', 'fontsize', font_size);


%% Save figure
% Prepare printer
h_fig.PaperPositionMode = 'auto';
h_fig.Units = 'Inches';
fig_pos = h_fig.Position;
set(h_fig, 'PaperUnits','Inches','PaperSize', [fig_pos(3), fig_pos(4)]);
% Set filename
output_filename = 'D_Grad_profile.pdf';
output_full_path = strcat(output_figures_folder, output_filename);
if bl_save_figures
    print(h_fig, output_full_path, '-dpdf', '-r0');
end;
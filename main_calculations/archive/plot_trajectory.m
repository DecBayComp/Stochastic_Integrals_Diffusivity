


fig_hand = figure(1);
set_my_fig_size (fig_hand);
clf;

plot(t_mesh, x_lambda, 'LineWidth', 2);
xlabel('t');
ylabel('x(t)');
title('Simulated trajectory for lambdas: 0(blue), 0.5(red) and 1(yellow)');

%% Saving
output_filename_png = sprintf('D%i_F%i_03_Simulated trajectories.png', D_case_number, f_case_number);
output_full_path_svg = strcat(output_figures_folder, output_filename_png);
% output_full_path_pdf = strcat(output_figures_folder, output_filename_pdf);
if bl_save_figures
    fig_hand.PaperPositionMode = 'auto';
    saveas(fig_hand, output_full_path_svg, 'png');
%     saveas(fig_hand, output_full_path_pdf, 'pdf');
    disp('Figure saved (trajectory)');
end;






















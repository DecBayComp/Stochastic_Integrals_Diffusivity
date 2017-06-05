function set_article_figure_size(h_fig, rows, cols, height_factor)
% function set_article_figure_size(h_fig, height_to_width_ratio, width_in_columns)

%% Globals
% % % global font_size;
% % % global one_column_figure_width_in;
% % % global mac_screen_dpi;

%% Constants
load_constants;
mac_screen_dpi = 120;
ubuntu_screen_dpi = 150;
height_base_in = 2.1;
screen_dpi = ubuntu_screen_dpi;
one_column_figure_width_in = 8.7/2.54;
two_columns_figure_width_in = 17.8/2.54;
one_third_page_figure_width_in = 17.8 / 2.54 / 3;
widths_array_in = [one_column_figure_width_in, two_columns_figure_width_in, one_third_page_figure_width_in];


%% Calculate
width = widths_array_in(cols);
height = height_base_in * rows * height_factor;


%% Adjust
paper_size = [width, height];
fig_size = paper_size * screen_dpi;
% Set screen size
h_fig.Units = 'pixels';
fig_pos = h_fig.Position;
fig_pos(3) = fig_size(1);
fig_pos(4) = fig_size(2);
h_fig.Position = fig_pos;
% Set font size
set(0, 'DefaultAxesFontSize', font_size);









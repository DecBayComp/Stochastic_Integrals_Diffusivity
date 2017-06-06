

function [selected_bins_indices, selected_bins_centers] = get_selected_bins_indices(data_struct)

%% Globals
% global x_bins_centers_saved;


%% Constants
load_constants;

%% Load data from the data structure
x_bins_centers = data_struct.x_bins_centers;

% Identifying bins numbers from the given coordinates
selected_bins_count = length(selected_x_over_L);
selected_bins_indices = zeros(1, selected_bins_count);
selected_bins_centers = zeros(1, selected_bins_count);
% % selected_bins_centers = zeros(1, selected_bins_count);
for i_bin = 1:selected_bins_count
    [~, min_number]= min((selected_x_over_L(i_bin) * L - x_bins_centers).^2);
    selected_bins_indices(i_bin) = min_number;
end;
selected_bins_centers(:) = x_bins_centers(selected_bins_indices);








function [selected_bins_indices, selected_bins_centers] = get_selected_bins_indices(data_struct)

%% Globals
% global x_bins_centers_saved;


%% Constants
load_constants;

%% Load data from the data structure
x_bins_centers_saved = data_struct.x_bins_centers_saved;

% Identifying bins numbers from the given coordinates
selected_bins_count = length(selected_x_over_L);
selected_bins_indices = zeros(lambda_count, selected_bins_count);
selected_bins_centers = zeros(lambda_count, selected_bins_count);
% % selected_bins_centers = zeros(1, selected_bins_count);
for l_ind = 1:lambda_count
    for i_bin = 1:selected_bins_count
        [~, min_number]= min((selected_x_over_L(i_bin) * L - x_bins_centers_saved{l_ind}).^2);
        selected_bins_indices(l_ind, i_bin) = min_number;
    end;
    selected_bins_centers(l_ind, :) = x_bins_centers_saved{l_ind}(selected_bins_indices(l_ind, :));
end;





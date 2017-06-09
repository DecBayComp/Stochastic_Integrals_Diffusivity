


function [x_bins_borders, x_bins_centers, bins_number, x_bins_widths,...
    points_count_in_bins, points_binned] = select_bins_adaptive_mesh(x_data, dx_data, points_in_bin)
%% The calculations assume no points have exactly the same locations


%% Constants
REL_PRECISION = 1e-4;


%% Initialize
N = length(x_data);
%= Calculate the number of bins
%= The excess of the points will go into the first bin, so I start binning
%= from the right-hand end
bins_number = floor(N/points_in_bin);
% Sort data
[x_data_sorted, sorted_indices] = sort(reshape(x_data, 1, []), 'ascend');
dx_data_sorted = reshape(dx_data(sorted_indices), 1, []);
% Determine zone size
x_min = x_data_sorted(1);
x_max = x_data_sorted(end);
% Initialize output containers
x_bins_borders = zeros(bins_number, 2);
points_count_in_bins = zeros(1, bins_number);
points_binned = cell(1, bins_number);


%% Calculate bin borders and bin the data
x_bins_borders(bins_number, 2) = x_max;
for bin = bins_number:-1:2
    % Find and save borders
    bin_from_end = bins_number - bin + 1;
    previous_border = (x_data_sorted(end - bin_from_end * points_in_bin + 1) + x_data_sorted(end - bin_from_end * points_in_bin)) / 2;
    x_bins_borders(bin, 1) = previous_border; 
    x_bins_borders(bin - 1, 2) = previous_border; 
    
    % Bin data
    indices = (N - bin_from_end * points_in_bin + 1) : (N - (bin_from_end - 1) * points_in_bin);
    points_binned{bin} = [x_data_sorted(indices); dx_data_sorted(indices)];
    points_count_in_bins(bin) = points_in_bin;
end
% The first bin includes all the remaining points
x_bins_borders(1, 1) = x_min;
indices = 1 : (N - (bins_number - 1) * points_in_bin);
points_binned{1} = [x_data_sorted(indices); dx_data_sorted(indices)];
points_count_in_bins(1) = length(indices);


%% Slightly shift the first and the last boundary to be sure to include the boundary points
x_bins_borders(1, 1) = x_bins_borders(1, 1) - (x_bins_borders(1, 2) - x_bins_borders(1, 1)) * REL_PRECISION;
x_bins_borders(end, 2) = x_bins_borders(end, 2) + (x_bins_borders(end, 2) - x_bins_borders(end, 1)) * REL_PRECISION;


%% Binning finished. Prepare the output
x_bins_widths = x_bins_borders(:, 2) - x_bins_borders(:, 1);
x_bins_centers = mean(x_bins_borders, 2);













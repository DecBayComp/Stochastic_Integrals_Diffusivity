%% Use an existing space tesselation to bin data with a given hard limit of points per bin


function [elements_in_bins_count, points_in_bins, bl_empty_bins] = bin_into_predefined_bins(x, dx, x_bins_borders, n_limit)

% Sort x
[~, sorted_indices] = sort(x, 'ascend');
sorted_data = [x(sorted_indices); dx(sorted_indices)];
x_bins_number = size(x_bins_borders, 2);


% Initialize output arrays
points_in_bins = cell(1, x_bins_number);
elements_in_bins_count = zeros(1, x_bins_number);
bl_empty_bins = zeros(1, x_bins_number);


% Bin data in pre-defined bins
left_border_point = 1;
for bin = 1:x_bins_number
    bl_empty_bin = false;

    % Get points falling into the current bin
    bl_in_bin = (sorted_data(1, :) > x_bins_borders(bin, 1)) & (sorted_data(1, :) <= x_bins_borders(bin, 2));
    jumps_indices = find(bl_in_bin);
    tot_points_in_bin = sum(bl_in_bin);

    % Impose a hard limit of points if necessary
    if n_limit > 0 && tot_points_in_bin > n_limit
        tot_points_in_bin_kept = n_limit;

        % Make a random draw from points indices
        sel_indices = randperm(tot_points_in_bin, n_limit);
        % Convert to jumps indices of the original array
        jumps_indices = jumps_indices(sel_indices);
    else
        tot_points_in_bin_kept = tot_points_in_bin;
    end

    % Bin points
    elements_in_bins_count(bin) = tot_points_in_bin_kept;
    right_border_point = left_border_point + tot_points_in_bin - 1;

    % Assign to bin
    points_in_bins{bin} = sorted_data(:, jumps_indices);

    % Move to next bin, left side
    left_border_point = right_border_point + 1;

    %% If bin empty for this trial, set flag
    if tot_points_in_bin_kept == 0
        bl_empty_bins(bin) = true;
    end

end
  
    
    





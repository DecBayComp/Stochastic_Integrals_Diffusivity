%% This function uses decimation (division in halfs) to create an adaptive mesh


function [x_bins_borders, x_bins_centers, x_bins_number, x_bins_widths,...
    elements_in_bins_count, binned_jumps] = select_bins_adaptive_mesh(x_data, dx_data, min_points_in_bin)

% %% Globals
% global N;


%% Constants
% load_constants;
% initial_bin_width_multiplier = 1/8;
% bin_width_increase_factor = 1.2;
% MAX_ITERATIONS = 100;
% max_bins_number = 1000;
REL_PRECISION = 1e-4;


%% Initialize
% bin_width_multiplier = initial_bin_width_multiplier;
% dx_abs_mean = mean(abs(dx_data));

% % % Remove the last point, because no jump is associated with it
% % % I am not sure why I do this...
% % x_data(end) = [];
N = length(x_data);
max_bins_number = floor(N/min_points_in_bin);
bl_bins_indivisible = ones(1, max_bins_number);
x_bins_borders = zeros(max_bins_number, 2);
points_count_in_bins = zeros(1, max_bins_number);
points_in_bins = cell(1, max_bins_number);


% Determine the size of the explored zone
x_min = min(x_data);
x_max = max(x_data);
zone_size = x_max - x_min;
% Slightly shifting the zone limits to avoid points at the border
x_max = x_max + zone_size * REL_PRECISION/2;
x_min = x_min - zone_size * REL_PRECISION/2;
% zone_size = x_max - x_min;


%% Initialize by setting the whole zone to the first bin
total_bins_number = 1;
x_bins_borders(1, :) = [x_min, x_max];
points_count_in_bins(1) = N;
bl_bins_indivisible(1) = 0;
% Sort the original points before assigning to the first bin
[~, sorted_indices] = sort(x_data, 'ascend');
points_in_bins{1} = [reshape(x_data(sorted_indices), 1, []); reshape(dx_data(sorted_indices), 1, [])];


%% Repeat the following operations until no divisions can be made
bl_flag_division_performed = 1;
while(bl_flag_division_performed)
    bl_flag_division_performed = 0;
    %% Attempt to divide the bin in two halves
    total_bins_number_snapshot = total_bins_number;
    for bin = 1:total_bins_number_snapshot
        % Skip bin if indivisible
        if bl_bins_indivisible(bin)
            continue;
        end;

        % Calculate borders of the new bins
        cur_width = x_bins_borders(bin, 2) - x_bins_borders(bin, 1);
        new_bin_1_borders = x_bins_borders(bin, 1) + [0, cur_width/2];
        new_bin_2_borders = new_bin_1_borders + cur_width/2;
        % Calculate the number of points in each bin
        points_count_in_new_bin_1 = sum(points_in_bins{bin}(1, :) <= new_bin_1_borders(2));
        points_count_in_new_bin_2 = points_count_in_bins(bin) - points_count_in_new_bin_1;
        % Accept the division if in both new bins the number of points is greater than the minimum
        if points_count_in_new_bin_1 < min_points_in_bin || points_count_in_new_bin_2 < min_points_in_bin
            bl_bins_indivisible(bin) = 1;
            continue;
        end;
        % Else distribute the points
        points_in_new_bin_1 = points_in_bins{bin}(:, 1:points_count_in_new_bin_1);
        points_in_new_bin_2 = points_in_bins{bin}(:, (points_count_in_new_bin_1 + 1):end);

        % Save bin_1
        points_in_bins{bin} = points_in_new_bin_1;
        points_count_in_bins(bin) = points_count_in_new_bin_1;
        x_bins_borders(bin, :) = new_bin_1_borders;
        bl_bins_indivisible(bin) = points_count_in_new_bin_1 < min_points_in_bin * 2;

        % Save bin_2
        total_bins_number = total_bins_number + 1;
        points_in_bins{total_bins_number} = points_in_new_bin_2;
        points_count_in_bins(total_bins_number) = points_count_in_new_bin_2;
        x_bins_borders(total_bins_number, :) = new_bin_2_borders;
        bl_bins_indivisible(total_bins_number) = points_count_in_new_bin_2 < min_points_in_bin * 2;

        % Flag that a division was performed
        bl_flag_division_performed = 1;
    end;
end;


%% Binning is completed, but I have to sort the bins
% Cut the bins that are not used
x_bins_borders((total_bins_number + 1):end, :) = [];
points_count_in_bins((total_bins_number + 1):end, :) = [];
points_in_bins((total_bins_number + 1):end) = [];


x_bins_centers = mean(x_bins_borders, 2);
[~, sorted_indices] = sort(x_bins_centers, 'ascend');
% Reorder results
x_bins_centers = x_bins_centers(sorted_indices);
x_bins_borders = x_bins_borders(sorted_indices, :);
points_count_in_bins = points_count_in_bins(sorted_indices);
points_in_bins = points_in_bins(sorted_indices);


%% Bins sorted. Prepare the output
x_bins_number = total_bins_number;
x_bins_widths = x_bins_borders(:, 2) - x_bins_borders(:, 1);
elements_in_bins_count = points_count_in_bins;    
binned_jumps = points_in_bins;

1;

    














function [x_bins_borders, x_bins_centers, x_bins_number, x_bins_widths,...
    elements_in_bins_count, binned_jumps] = select_bins_equal_width(x_data, dx_data, min_points_in_bin)

%% Globals
% global min_points_in_bin;


%% Constants
initial_bin_width_multiplier = 1/8;
bin_width_increase_factor = 1.2;
MAX_ITERATIONS = 100;
max_bins_number = 1000;
REL_PRECISION = 1e-4;


%% Initialize
bin_width_multiplier = initial_bin_width_multiplier;
dx_abs_mean = mean(abs(dx_data));
% Remove the last point, because it's only an end. No jump associated
x_data(end) = [];


% Determine the size of the explored zone
x_min = min(x_data);
x_max = max(x_data);
zone_size = x_max - x_min;
% Slightly shifting the zone limits to avoid points at the border
x_max = x_max + zone_size * REL_PRECISION/2;
x_min = x_min - zone_size * REL_PRECISION/2;
zone_size = x_max - x_min;
    

%% Try the smallest bin size first and increase it if it doesn't fit the min. points requirement
for try_number = 0:MAX_ITERATIONS-1
    %% Setting the bin size
    x_bins_number = round(zone_size / (dx_abs_mean * bin_width_multiplier));
    if x_bins_number > max_bins_number
        x_bins_number = max_bins_number;
    end;
    x_bin_width = zone_size / x_bins_number;
    x_bins_borders = x_min + (0:x_bins_number) * x_bin_width;
    x_bins_centers = (x_bins_borders(1:end-1) + x_bins_borders(2:end))/2;


    %% Binning the trajectory 
    % Calculating the number of elements
    bin_occupation_numbers_raw = ceil((x_data - x_min) / x_bin_width);
    elements_in_bins_count = histcounts(bin_occupation_numbers_raw, 0.5 + (0:x_bins_number));   
    % 0.5 makes sure the numbers fall within bins and not on the borders (and in the right bins)

    
    %% Checking for empty bins 
    if sum(elements_in_bins_count < min_points_in_bin) == 0
        % Only if all bins have enough points, continue the calculations
        break;
    end;
    % Otherwise adjust the bin size and repeat
    bin_width_multiplier = bin_width_multiplier * bin_width_increase_factor;

end;


%% Bins are selected. Bin the trajectory
binned_jumps = cell(1, x_bins_number);
% Each cell contains elements in the form (x, dx)
for bin = 1:x_bins_number
    indices = find(bin_occupation_numbers_raw == bin);
%     disp(indices);
%     fprintf('\n');
    binned_jumps{bin} = [x_data(indices); dx_data(indices)];
end;


x_bins_widths = x_bin_width * ones(x_bins_number, 1);












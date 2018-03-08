%% Create an adaptive mesh of the space.
% Conditions:
% 1. Approximately same number of points per bin
% 2. The minimum bin size is related to the mean jump size in this bin through min_bin_to_jump_ratio
%
% The max_points_per_bin will impose a hard limit on the number of points by randomly choosing after binning. 
% The value can be a vector, then several configurations will be generated.
% -1 means no limit



function [x_bins_borders, x_bins_centers, bins_number, x_bins_widths,...
    point_count_in_bins, variance_in_bins, points_binned] = select_bins_adaptive_mesh(x_data, dx_data, max_points_per_bin)
%% The calculations assume no points have exactly the same locations


% Start timer
tic;

%% Constants
load_constants;
REL_PRECISION = 1e-4;
initial_bins_number = 100;
N_for_binning = 1e6;


conv_factor = 1/2;	% convergence rate to the desired bin width


%% Initialize
N = length(x_data);

% To save memory on large datasets, drop some of the points to optimize performance
if N > N_for_binning
    % Select indices
    sel_indices = randperm(N, N_for_binning);
    
    % Reduce data size
    x_data = x_data(sel_indices);
    dx_data = dx_data(sel_indices);
    N = N_for_binning;
end


% The excess of the points will go into the first bin, so start binning
% from the right-hand end
bins_number = initial_bins_number;
points_per_bin = floor(N / bins_number);

% Sort data
[x_data_sorted, sorted_indices] = sort(reshape(x_data, 1, []), 'ascend');
dx_data_sorted = reshape(dx_data(sorted_indices), 1, []);

% Determine zone size
x_min = x_data_sorted(1);
x_max = x_data_sorted(end);

% Initialized bins borders in x and point number
x_bins_borders = zeros(bins_number, 2);
indices_bins_borders = zeros(bins_number, 2);

% Initialize output containers
point_count_in_bins = zeros(1, bins_number);
points_binned = cell(1, bins_number);
variance_in_bins = zeros(1, bins_number);
mean_jumps = zeros(1, bins_number);


%% Calculate bin borders and bin the data
% Set the right border of the last bin
% Indices borders are always included into their bin (the indexes of points
% at the borders of the bins)
x_bins_borders(bins_number, 2) = x_max;
indices_bins_borders(bins_number, 2) = N;
first_bin = 1;
for bin = bins_number:-1:1
	% Print progress
	fprintf('Binning: Processing bin %i/%i.\n', bins_number - bin + 1, bins_number);
	
	% Store right border index (fixed)
	right_border_ind = indices_bins_borders(bin, 2);
	
    % Estimate the left bin border based on the point number
	left_border_ind = right_border_ind - points_per_bin + 1;
	
	% Keep increasing the bin width, while the mean jump is smaller than a
	% pre-defined bin width fraction and while there are unbinned points left
	bl_small_mean_jump = false;
	ind = left_border_ind;
    while ind >= 1  % try to find the index of the jump which would make the left border
		% Calculate the mean jump in bin
		mean_jump_length = std(dx_data_sorted(ind:right_border_ind));
		
		% Calculate bin width
		bin_width = x_data_sorted(right_border_ind) - x_data_sorted(ind);
		
		% Save and leave cycle if the condition is satisfied
        if bin_width >= mean_jump_length * min_bin_to_jump_ratio
			% Save borders
			indices_bins_borders(bin, 1) = ind;
            if ind == 1
                % If last point
                x_bins_borders(bin, 1) = x_data_sorted(ind);
            else
                x_bins_borders(bin, 1) = (x_data_sorted(ind) + x_data_sorted(ind - 1)) / 2;
            end
            
            % Change the bin to the left if not the last one
            if bin > 1
                indices_bins_borders(bin - 1, 2) = ind - 1;
                x_bins_borders(bin - 1, 2) = x_bins_borders(bin, 1);
            end
			
			% Bin the data
			indices = ind:right_border_ind;
			points_binned{bin} = [x_data_sorted(indices); dx_data_sorted(indices)];
			point_count_in_bins(bin) = right_border_ind - ind + 1;
			
			% Calculate varaince
			variance_in_bins(bin) = var(points_binned{bin}(2, :));
			mean_jumps(bin) = mean_jump_length;
			
			% Print progress
			fprintf('The bin was smaller than %.1f * (mean jump) and was enlarged by %i points.\n', min_bin_to_jump_ratio, right_border_ind - ind + 1 - points_per_bin);
			
			% Set flag
			bl_small_mean_jump = true;
			break;
        end
		
		%% Estimate the index of the border which would give the right bin size
		est_width_increase_ratio = mean_jump_length * min_bin_to_jump_ratio / bin_width;
		
		% Estimated number of points to add
		est_points_increase = (right_border_ind - ind + 1) * (est_width_increase_ratio - 1);
		
		% Reduce by a convergence factor to be more accurate
		est_points_increase = ceil(est_points_increase * conv_factor) + 1;
		
		% Calculate the next index to consider
		ind = ind - est_points_increase;
		1;
		
    end
	
	% If the condition was not satisfied, then we ran out of points, and the current points will be attached to the previous bin
	if ~bl_small_mean_jump 
		indices_bins_borders(bin + 1, 1) = 1;
		x_bins_borders(bin + 1, 1) = x_data_sorted(1);
		
		% Rebin the data in the previous bin
		right_border_ind = indices_bins_borders(bin + 1, 2);
		left_border_ind = 1;
		indices = left_border_ind:right_border_ind;
		points_binned{bin + 1} = [x_data_sorted(indices); dx_data_sorted(indices)];
		point_count_in_bins(bin + 1) = right_border_ind - left_border_ind + 1;

		% Recalculate varaince
		variance_in_bins(bin + 1) = var(points_binned{bin + 1}(2, :));
		
		% Reduce total bins number
		first_bin = bin + 1;
		bins_number = bins_number - first_bin + 1;
		
		% Print progress
		fprintf('Binning: Ran out of points. Setting up the last bin #%i\n', bins_number);
		
		% Exit bin cycle
		break;
	end
end



%% Remove empty bins in the start (if any)
x_bins_borders = x_bins_borders(first_bin:end, :);
indices_bins_borders = indices_bins_borders(first_bin:end, :);
variance_in_bins = variance_in_bins(first_bin:end);
mean_jumps = mean_jumps(first_bin:end);
point_count_in_bins = point_count_in_bins(first_bin:end);

% Binned points
points_binned_new = cell(1, bins_number);
for bin = 1:bins_number
	points_binned_new{bin} = points_binned{first_bin + bin - 1};
end
points_binned = points_binned_new;



%% Print some statistics
fprintf('Binning: %i bins were smaller than %.1f * (mean jump) and were enlarged.\n', sum(point_count_in_bins > points_per_bin),min_bin_to_jump_ratio);



%% Slightly shift the first and the last boundary to be sure to include the boundary points
x_bins_borders(1, 1) = x_bins_borders(1, 1) - (x_bins_borders(1, 2) - x_bins_borders(1, 1)) * REL_PRECISION;
x_bins_borders(end, 2) = x_bins_borders(end, 2) + (x_bins_borders(end, 2) - x_bins_borders(end, 1)) * REL_PRECISION;


%% Binning finished. Impose the hard limit of points
n_limits = max_points_per_bin;
n_limits_count = length(max_points_per_bin);
points_binned_new = cell(n_limits_count, bins_number);
variance_in_bins_new = zeros(n_limits_count, bins_number);
point_count_in_bins_new = zeros(n_limits_count, bins_number);

for lim_ind = 1:n_limits_count
    n_limit = n_limits(lim_ind);
    
    % For each bin create a version with different number of points per bin
    for bin = 1:bins_number
        tot_points_in_bin = point_count_in_bins(bin);

        % Create a random draw from point indices in bin (if there are more points than required)
        if n_limit > 0 && n_limit < tot_points_in_bin
            sel_indices = randperm(tot_points_in_bin, n_limit);
        else
            sel_indices = 1:tot_points_in_bin;
        end
        
        % Store the selected points
        points_binned_new{lim_ind, bin} = points_binned{bin}(:, sel_indices);
        
        % Calculate variance
        variance_in_bins_new(lim_ind, bin) = var(points_binned_new{lim_ind, bin}(2, :));
        
        % Count points
        point_count_in_bins_new(lim_ind, bin) = length(sel_indices);
    end
end

% Rewrite variables to store and return results
points_binned = points_binned_new;
variance_in_bins = variance_in_bins_new;
point_count_in_bins = point_count_in_bins_new;



%% Binning finished. Prepare the output
x_bins_widths = x_bins_borders(:, 2) - x_bins_borders(:, 1);
x_bins_centers = mean(x_bins_borders, 2);

% Print execution time
fprintf('Binning: Completed in %.2f min\n', toc/60);

1;











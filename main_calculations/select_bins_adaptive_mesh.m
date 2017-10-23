%% Create an adaptive mesh of the space.
% Conditions:
% 1. Approximately same number of points per bin
% 2. The minimum bin size is related to the mean jump size in this bin through min_bin_to_jump_ratio



function [x_bins_borders, x_bins_centers, bins_number, x_bins_widths,...
    points_count_in_bins, variance_in_bins, points_binned] = select_bins_adaptive_mesh(x_data, dx_data, points_in_bin)
%% The calculations assume no points have exactly the same locations


% Start timer
tic;

%% Constants
load_constants;
REL_PRECISION = 1e-4;

% The bin size will be increased by the following number of points each time that the bin is smaller than the mean jump
increase_step = ceil(points_in_bin / 1000);

conv_factor = 1/2;	% convergence rate to the desired bin width


%% Initialize
N = length(x_data);
% Calculate the number of bins
% The excess of the points will go into the first bin, so I start binning
% from the right-hand end
bins_number = floor(N/points_in_bin);

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
points_count_in_bins = zeros(1, bins_number);
points_binned = cell(1, bins_number);
variance_in_bins = zeros(1, bins_number);
mean_jumps = zeros(1, bins_number);


%% Calculate bin borders and bin the data
% Set the right border of the last bin
% Indices borders are always included into their bin
x_bins_borders(bins_number, 2) = x_max;
indices_bins_borders(bins_number, 2) = N;
for bin = bins_number:-1:1
	% Print progress
	fprintf('Binning: Processing bin %i/%i.\n', bins_number - bin + 1, bins_number);
	
	% Store rigth border index (fixed)
	right_border_ind = indices_bins_borders(bin, 2);
	
    % Estimate the left bin border based on the point number
	left_border_ind = right_border_ind - points_in_bin + 1;
	
	% Keep increasing the bin width, while the mean jump in bin is not small enough compared to the bins and while there are unbinned points left
	bl_small_mean_jump = false;
	ind = left_border_ind;
	while ind >= 1
		% Calculate the mean jump in bin
		mean_jump_length = std(dx_data_sorted(ind:right_border_ind));
		
		% Calculate bin width
		bin_width = x_data_sorted(right_border_ind) - x_data_sorted(ind);
		
		% Leave cycle if the condition is satisfied
		if bin_width >= mean_jump_length * min_bin_to_jump_ratio
			% Save borders
			indices_bins_borders(bin, 1) = ind;
			indices_bins_borders(bin - 1, 2) = ind - 1;
			x_bins_borders(bin, 1) = (x_data_sorted(ind) + x_data_sorted(ind - 1)) / 2;
			x_bins_borders(bin - 1, 2) = x_bins_borders(bin, 1);
			
			% Bin the data
			indices = ind:right_border_ind;
			points_binned{bin} = [x_data_sorted(indices); dx_data_sorted(indices)];
			points_count_in_bins(bin) = right_border_ind - ind + 1;
			
			% Calculate varaince
			variance_in_bins(bin) = var(points_binned{bin}(2, :));
			mean_jumps(bin) = mean_jump_length;
			
			% Print progress
			fprintf('The bin was smaller than %.1f * (mean jump) and was enlarged by %i points.\n', min_bin_to_jump_ratio, right_border_ind - ind + 1 - points_in_bin);
			
			% Set flag
			bl_small_mean_jump = true;
			break;
		end;
		
		%% Estimate the index of the border which would give the right bin size
		est_width_increase_ratio = mean_jump_length * min_bin_to_jump_ratio / bin_width;
		
		% Estimated number of points to add
		est_points_increase = (right_border_ind - ind + 1) * (est_width_increase_ratio - 1);
		
		% Reduce by a convergence factor to be more accurate
		est_points_increase = ceil(est_points_increase * conv_factor) + 1;
		
		% Calculate the next index to consider
		ind = ind - est_points_increase;
		
	end;
	
	% If the condition was not satisfied, then we ran out of points, and the current points will be attached to the previous bin
	if ~bl_small_mean_jump 
		indices_bins_borders(bin + 1, 1) = 1;
		x_bins_borders(bin + 1, 1) = x_data_sorted(1);
		
		% Rebin the data in the previous bin
		right_border_ind = indices_bins_borders(bin + 1, 2);
		left_border_ind = 1;
		indices = left_border_ind:right_border_ind;
		points_binned{bin + 1} = [x_data_sorted(indices); dx_data_sorted(indices)];
		points_count_in_bins(bin + 1) = right_border_ind - left_border_ind + 1;

		% Recalculate varaince
		variance_in_bins(bin + 1) = var(points_binned{bin + 1}(2, :));
		
		% Reduce total bins number
		first_bin = bin + 1;
		bins_number = bins_number - first_bin + 1;
		
		% Print progress
		fprintf('Binning: Ran out of points. Setting up the last bin #%i\n', bins_number);
		
		% Exit bin cycle
		break;
	end;
end



%% Remove empty bins in the start
x_bins_borders = x_bins_borders(first_bin:end, :);
indices_bins_borders = indices_bins_borders(first_bin:end, :);
variance_in_bins = variance_in_bins(first_bin:end);
mean_jumps = mean_jumps(first_bin:end);
points_count_in_bins = points_count_in_bins(first_bin:end);

% Binned points
points_binned_new = cell(1, bins_number);
for bin = 1:bins_number
	points_binned_new{bin} = points_binned{first_bin + bin - 1};
end;
points_binned = points_binned_new;



%% Print some statistics
fprintf('Binning: %i bins were smaller than %.1f * (mean jump) and were enlarged.\n', sum(points_count_in_bins > points_in_bin),min_bin_to_jump_ratio);



%% Slightly shift the first and the last boundary to be sure to include the boundary points
x_bins_borders(1, 1) = x_bins_borders(1, 1) - (x_bins_borders(1, 2) - x_bins_borders(1, 1)) * REL_PRECISION;
x_bins_borders(end, 2) = x_bins_borders(end, 2) + (x_bins_borders(end, 2) - x_bins_borders(end, 1)) * REL_PRECISION;


%% Binning finished. Prepare the output
x_bins_widths = x_bins_borders(:, 2) - x_bins_borders(:, 1);
x_bins_centers = mean(x_bins_borders, 2);

% Print execution time
fprintf('Binning: Completed in %.2f min\n', toc/60);

1;











% Use an existing space tesselation to bin data


function [] = bin_predefined_bins()


%% Bin points for the current trial
    x = trials_x(:, trial)';
    dx = trials_dx(:, trial)';
    cur_dx_bck_data = zeros(size(dx));    % deprecated
    % Sort x data
    [~, sorted_indices] = sort(x, 'ascend');
    sorted_data = [x(sorted_indices); dx(sorted_indices); cur_dx_bck_data(sorted_indices)];
    % Initialize binned data cell array
    points_in_bins = cell(1, x_bins_number);
    elements_in_bins_count = zeros(1, x_bins_number);
    % Bin data in pre-defined bins
    counter_start = 1;
	for bin = 1:x_bins_number
		% Initialize
		bl_empty_bin = false;
		
        elements_in_bins_count(bin) = sum((sorted_data(1, :) > x_bins_borders(bin, 1)) & (sorted_data(1, :) <= x_bins_borders(bin, 2)));
        counter_end = counter_start + elements_in_bins_count(bin) - 1;
        points_in_bins{bin} = sorted_data(:, counter_start : counter_end);
        counter_start = counter_end + 1;

% % %         %% If requested, keep only the minimum number of points per bin
% % % 		if bl_keep_only_min_points_in_bin && elements_in_bins_count(bin) > 
% % % 			cur_points = points_in_bins{bin};
% % % 			cur_points_count = elements_in_bins_count(bin);
% % % 			
% % % 		end
		
		%% Save calculated values to the data structure
        data_struct.n_j(bin) = elements_in_bins_count(bin);
		data_struct.dx_mean_in_bins(bin) = mean(points_in_bins{bin}(2, :));
        data_struct.V_j(bin) = var(points_in_bins{bin}(2, :));
		data_struct.mean_jump_length_bins = sqrt(data_struct.V_j(bin));
		data_struct.dx_mean_all_bins = mean(sorted_data(2, :));
		data_struct.V_all_bins = var(sorted_data(2, :));





%% Calculate Kolmogorov-Smirnov distance for a and b for each trial


function [trials_a_KS_distance, trials_b_KS_distance] = batch_calculate_KS_distance(trials, x_bins_number, data_struct, trials_data, trials_MAP_a, trials_MAP_b)

%% Constants
load_constants;
t_step = t_step;	% Make it visible to the parallel pool



%% Pre-calculate MAP distributions for each bin, simulation type and inference type
%% Initialize
% Trials indices for different simulation types
trials_indices = cell(lambda_types_count, 1);

% MAP distributions
MAP_a_distr_bin_sorted = cell(lambda_types_count, x_bins_number, conventions_count);
MAP_b_distr_bin_sorted = cell(lambda_types_count, x_bins_number);

% Calculate
for lambda_type = 1:lambda_types_count
	% Get all indices of trials for this simulation type
	trials_indices{lambda_type} = find(data_struct.trial_simulation_type == lambda_type);
	
	% Extract MAP a for all bins and conventions for this simulation type
	MAP_a_distr_bin = trials_MAP_a(trials_indices{lambda_type}, :, :, 1);
	
	% Drop NaN and sort MAP a
	for bin = 1:x_bins_number
		for convention = 1:conventions_count
			MAP_a_distr_bin_sorted{lambda_type, bin, convention} = MAP_a_distr_bin(:, bin, convention);
			MAP_a_distr_bin_sorted{lambda_type, bin, convention} = MAP_a_distr_bin_sorted{lambda_type, bin, convention}(~isnan(MAP_a_distr_bin_sorted{lambda_type, bin, convention}));
			MAP_a_distr_bin_sorted{lambda_type, bin, convention} = sort(MAP_a_distr_bin_sorted{lambda_type, bin, convention});
		end;
	end;
	
	% Extract MAP a for all bins for this simulation type
	MAP_b_distr_bin = trials_MAP_b(trials_indices{lambda_type}, :, 1);
	
	% Drop NaN and sort MAP b
	for bin = 1:x_bins_number
		MAP_b_distr_bin_sorted{lambda_type, bin} = MAP_b_distr_bin(:, bin);
		MAP_b_distr_bin_sorted{lambda_type, bin} = MAP_b_distr_bin_sorted{lambda_type, bin}(~isnan(MAP_b_distr_bin_sorted{lambda_type, bin}));
		MAP_b_distr_bin_sorted{lambda_type, bin} = sort(MAP_b_distr_bin_sorted{lambda_type, bin});
	end;
end;
	


% Initialize distances
trials_b_KS_distance = zeros(trials, x_bins_number);
trials_a_KS_distance = zeros(trials, x_bins_number, conventions_count);



%% Calculate distance
parfor trial = 1:trials
	fprintf('Calculating Kolmogorov-Smirnov distance. Trial: %i/%i\n', trial, trials);	
	
	% Identify simulation type
	lambda_type = data_struct.trial_simulation_type(trial);
	
	for bin = 1:x_bins_number
		
		% Load trial data
		cur_data_struct = trials_data{trial};
		
		% Skip empty bins
		if cur_data_struct.bl_empty_bin(bin)
			trials_b_KS_distance(trial, bin) = NaN;
			trials_a_KS_distance(trial, bin, :) = NaN;
			continue;
		end;
		
		
		
		%% Diffusivity b

		% Load b posterior for this bin and trial
		b_posterior_func_wrap = @(b) bin_b_posterior_func (bin, b, t_step, cur_data_struct, 'forward');

		% Calculate KS distance
		trials_b_KS_distance(trial, bin) = calculate_KS_distance(MAP_b_distr_bin_sorted{lambda_type, bin}, b_posterior_func_wrap);
		
		
		
	end;
end;

1;




		% Create temporary variable for loop parallelization
% 		tmp_trials_b_KS_distance = zeros(

			

			

% 			%% Drift a
% 
%  			%% Ito
% 			convention = 1;
% 			% Load the posterior for this bin and trial
% 			a_posterior_func_wrap = @(a) bin_a_posterior_func (bin, a, t_step, cur_data_struct, 'forward');
% 
% 			% Calculate distance
% 			trials_a_KS_distance(trial, convention, bin) = calculate_KS_distance(MAP_a_distr_bin, a_posterior_func_wrap);

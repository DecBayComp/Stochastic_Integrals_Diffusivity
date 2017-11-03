%% Calculate Kolmogorov-Smirnov distance for a and b for each trial


function [trials_a_KS_distance, trials_b_KS_distance] = batch_calculate_KS_distance(trials, trials_MAP_a, trials_MAP_b)

%% Constants
load_constants;



%% Pre-calculate the MAP distributions for each bin, simulation type and inference type







% MAP b distribution over trials and b posterior from one trial
trials_b_KS_distance = zeros(trials, x_bins_number);
trials_a_KS_distance = zeros(trials, x_bins_number, conventions_count);
trials_indices = cell(lambda_types_count, 1);

for lambda_type = 1:lambda_types_count
	%% Initialize
	
	% Get all indices of different simulation types
	trials_indices{lambda_type} = find(data_struct.trial_simulation_type == lambda_type);
	
	% Extract MAP a for all bins and conventions for this simulation type
	MAP_a_distr_bin = trials_MAP_a(trials_indices{lambda_type}, :, :, 1);
	
	% Drop NaN and sort MAP a
	MAP_a_distr_bin_sorted = cell(x_bins_number, conventions_count);
	for bin = 1:x_bins_number
		for convention = 1:conventions_count
			MAP_a_distr_bin_sorted{bin, convention} = MAP_a_distr_bin(:, bin, convention);
			MAP_a_distr_bin_sorted{bin, convention} = MAP_a_distr_bin_sorted{bin, convention}(~isnan(MAP_a_distr_bin_sorted{bin, convention}));
			MAP_a_distr_bin_sorted{bin, convention} = sort(MAP_a_distr_bin_sorted{bin, convention});
		end;
	end;
	
	% Extract MAP a for all bins for this simulation type
	MAP_b_distr_bin = trials_MAP_b(trials_indices{lambda_type}, :, 1);
	
	% Drop NaN and sort MAP b
	MAP_b_distr_bin_sorted = cell(x_bins_number, 1);
	for bin = 1:x_bins_number
		MAP_b_distr_bin_sorted{bin} = MAP_b_distr_bin(:, bin);
		MAP_b_distr_bin_sorted{bin} = MAP_b_distr_bin_sorted{bin}(~isnan(MAP_b_distr_bin_sorted{bin}));
		MAP_b_distr_bin_sorted{bin} = sort(MAP_b_distr_bin_sorted{bin});
	end;
	
	%% Start processing for the current simulation type
	parfor bin = 1:x_bins_number
		% Create temporary variable for loop parallelization
% 		tmp_trials_b_KS_distance = zeros(
		for trial = trials_indices{lambda_type}
			fprintf('Calculating Kolmogorov-Smirnov distance for %s simulations. Bin: %i/%i. Trial: %i/%i\n', lambda_types_names{lambda_type}, bin, x_bins_number, trial, trials);	
			% Load trial data
			cur_data_struct = trials_data{trial};

			% Skip empty bins
			if cur_data_struct.bl_empty_bin(bin)
				trials_b_KS_distance(trial, bin) = NaN;
% 				trials_a_KS_distance(trial, bin, :) = NaN;
				continue;
			end;

			%% Diffusivity b

			% Load the posterior for this bin and trial
			b_posterior_func_wrap = @(b) bin_b_posterior_func (bin, b, t_step, cur_data_struct, 'forward');

			% Calculate distance
			trials_b_KS_distance(trial, bin) = calculate_KS_distance(MAP_b_distr_bin_sorted{bin}, b_posterior_func_wrap);

% 			%% Drift a
% 
%  			%% Ito
% 			convention = 1;
% 			% Load the posterior for this bin and trial
% 			a_posterior_func_wrap = @(a) bin_a_posterior_func (bin, a, t_step, cur_data_struct, 'forward');
% 
% 			% Calculate distance
% 			trials_a_KS_distance(trial, convention, bin) = calculate_KS_distance(MAP_a_distr_bin, a_posterior_func_wrap);

		end;
	end;
end;
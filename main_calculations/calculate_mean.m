%% Calculate mean for each simulation type and n limit separately



function stat_struct = calculate_mean(stat_struct)



MAP_D_mean = zeros(ksi_count, x_bins_number, 4);
MAP_b_mean = zeros(ksi_count, x_bins_number, 4);
MAP_bb_prime_regular_interp_mean = zeros(ksi_count, x_bins_number);
n_j_mean = zeros(ksi_count, x_bins_number);
log_K_L_mean = zeros(ksi_count, conventions_count);
log_K_L_std = zeros(ksi_count, conventions_count);

for ksi_ind = 1:ksi_count
    % Mean
    MAP_D_mean(ksi_ind, :, :) = mean(trials_MAP_D(trials_ksi_type == ksi_ind, :, :), 1, 'omitnan' );
	MAP_b_mean(ksi_ind, :, :) = mean(trials_MAP_b(trials_ksi_type == ksi_ind, :, :), 1, 'omitnan' );
%     MAP_a_mean(lambda_type, :, :, :) = mean(trials_MAP_a(trial_simulation_type == lambda_type, :, :, :), 1, 'omitnan' );
    MAP_bb_prime_regular_interp_mean(ksi_ind, :) = mean(trials_MAP_bb_prime_regular_interp(trials_ksi_type == ksi_ind, :), 1, 'omitnan' );
	
	% n_j
	n_j_mean(ksi_ind, :, :) = mean(trials_n_j(trials_ksi_type == ksi_ind, :), 1, 'omitnan' );
        
    for convention = 1:conventions_count
        vals = reshape(trials_log_K_L(trials_ksi_type == ksi_ind, middle_bin, convention),[], 1);
        log_K_L_mean(ksi_ind, convention) = mean(vals, 'omitnan');
        log_K_L_std(ksi_ind, convention) = std(vals, 'omitnan');
    end
end

% Save
data_struct.MAP_D_mean = MAP_D_mean;
data_struct.MAP_b_mean = MAP_b_mean;
data_struct.MAP_bb_prime_regular_interp_mean = MAP_bb_prime_regular_interp_mean;
data_struct.n_j_mean = n_j_mean;
data_struct.mean_log_K_L = log_K_L_mean;
data_struct.std_log_K_L = log_K_L_std;
data_struct.eb_log_K_L = log_K_L_std * sqrt(2) * erfinv(0.95);























end
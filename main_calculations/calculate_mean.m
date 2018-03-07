%% Calculate mean for each simulation type and n limit separately



function stat_struct = calculate_mean(stat_struct)

% Initialize
ksi_count = length(stat_struct.ksi_array);
trials_ksi_type = stat_struct.trials_ksi_type;
n_limits_count = size(stat_struct.MAP_b, 2);
x_bins_number = size(stat_struct.MAP_b, 3);
conventions_count = size(stat_struct.log_K_L, 4);


MAP_D_mean = zeros(ksi_count, n_limits_count, x_bins_number, 4);
MAP_b_mean = zeros(ksi_count, n_limits_count, x_bins_number, 4);
MAP_bb_prime_regular_interp_mean = zeros(ksi_count, n_limits_count, x_bins_number);
n_j_mean = zeros(ksi_count, n_limits_count, x_bins_number);
log_K_L_mean = zeros(ksi_count, n_limits_count, conventions_count);
log_K_L_std = zeros(ksi_count, n_limits_count, conventions_count);

for ksi_ind = 1:ksi_count
    for lim_ind = 1:n_limits_count
        % Mean
        MAP_D_mean(ksi_ind, lim_ind, :, :) = mean(stat_struct.MAP_D(trials_ksi_type == ksi_ind, lim_ind, :, :), 1, 'omitnan' );
        MAP_b_mean(ksi_ind, lim_ind, :, :) = mean(stat_struct.MAP_b(trials_ksi_type == ksi_ind, lim_ind, :, :), 1, 'omitnan' );
        MAP_bb_prime_regular_interp_mean(ksi_ind, n_limits_count, :) = mean(stat_struct.MAP_bb_prime_regular_interp(trials_ksi_type == ksi_ind, lim_ind, :), 1, 'omitnan' );

        % n_j
        n_j_mean(ksi_ind, lim_ind, :, :) = mean(stat_struct.n_j(trials_ksi_type == ksi_ind, lim_ind, :), 1, 'omitnan' );

        for convention = 1:conventions_count
            vals = reshape(stat_struct.log_K_L(trials_ksi_type == ksi_ind, lim_ind, stat_struct.middle_bin, convention),[], 1);
            log_K_L_mean(ksi_ind, lim_ind, convention) = mean(vals, 'omitnan');
            log_K_L_std(ksi_ind, lim_ind, convention) = std(vals, 'omitnan');
        end
    end
end

% Save
stat_struct.MAP_D_mean = MAP_D_mean;
stat_struct.MAP_b_mean = MAP_b_mean;
stat_struct.MAP_bb_prime_regular_interp_mean = MAP_bb_prime_regular_interp_mean;
stat_struct.n_j_mean = n_j_mean;
stat_struct.log_K_L_mean = log_K_L_mean;
stat_struct.log_K_L_std = log_K_L_std;
stat_struct.log_K_L_eb = log_K_L_std * sqrt(2) * erfinv(0.95);























end
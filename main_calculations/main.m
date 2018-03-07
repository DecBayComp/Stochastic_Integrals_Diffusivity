set (0, 'DefaultAxesFontSize', 12);


%% Clear the workspace (do not use if want to be able to use bl_reload_trajectories = false)
% clear;

%% Constants
bl_force = true;
load_constants;

KSI_PRECISION = 1e-2;
bl_reload_trajectories = true;


%% Initialize
low_conf_percentile = (1 - CONF_LEVEL)/2;
high_conf_percentile = 1 - low_conf_percentile;
fD_min = [];
fD_max = [];
fD_step = [];
fD_mesh = [];
fD_pdf_plot_data = [];
pdf_norm = [];



%% Identify if calculation with force or with no force is required
% % % if bl_force
% % % % 	input_data_folder = '/home/aserov/Documents/Calculated_data/dilemma_with_force/';
% % % 	selected_f_case = enum_force_case;
% % % 	str_force = 'with_force';	
% % % else
% % % % 	input_data_folder = '/home/aserov/Documents/Calculated_data/dilemma_no_force/';
% % % 	selected_f_case = enum_no_force_case;
% % % 	str_force = 'no_force';
% % % end


%% Access trajectories folder and start loading trajectories
% Count the number of csv trajectories in a folder
cur_dir = dir([input_data_folder, '*.csv']);
% input_files_count = sum(~[cur_dir.isdir]);
input_files_count = 31;



% Load trajectories
if bl_reload_trajectories
    % Initialize arrays
    trials_x = zeros(N, input_files_count);
    trials_dx = zeros(N, input_files_count);
%     trials_lambdas = zeros(1, input_files_count);
	trials_D_case = zeros(1, input_files_count);
    trials_ksi = zeros(1, input_files_count);
% 	trials_f_case = zeros(1, input_files_count);
    % Load
    parfor file_num = 1:input_files_count
        filename = sprintf('sim_data_%09i.csv', file_num);
        fprintf('Loading trajectory from  file %s. Progress: %i/%i\n', filename, file_num, input_files_count);
        output_full_path = strcat(input_data_folder, filename);

        input_data = dlmread(output_full_path, CSV_DELIMITER);
		
		trials_D_case(file_num) = input_data(1,1);
		trials_ksi(file_num) = round(input_data(1,2), -log10(KSI_PRECISION));
        trials_x(:, file_num) = input_data(2:N+1, 1);
        trials_dx(:, file_num) = input_data(2:N+1,2);
    end
    input_data = [];
end



% %% Only keep data that corresponds to the currently analyzed condition (force / no force)
% trials_x = trials_x(:, trials_f_case == selected_f_case);
% trials_dx = trials_dx(:, trials_f_case == selected_f_case);
trials = input_files_count;
trials_data = cell(1, trials);

n_limits_count = length(n_limits);
	


% Identify and enumerate simulated ksi values
[ksi_array, ksi_count, trials_ksi_type, trial_first_ksi_type_index] = identify_ksi(trials_ksi);


%% Identify suitable bin locations based on all points for all trials
[x_bins_borders, x_bins_centers, x_bins_number, x_bins_widths,...
            ~, variance_in_bins, ~] = select_bins_adaptive_mesh(trials_x(:), trials_dx(:), -1);

% Estimate dx_Mean and V used for prior only. Average over everything
dx_mean_all_bins_all_trials = mean(trials_dx(:));
V = var(trials_dx(:));

% Calculate mean jump in each bin over all trials
mean_jump_bins_all_trials = sqrt(variance_in_bins);


%% Prepare the fine mesh
x_theor_min = x_min;
x_theor_max = x_max;
x_theor_step = (x_theor_max - x_theor_min) / (fine_mesh_steps_count - 1);
x_fine_mesh = x_theor_min:x_theor_step:x_theor_max;


%% Calculate true D and fD values on a fine mesh and in bins
% Calculate on the fine mesh
[D_theor_fine_data, D_grad_theor_fine_data] = D_func(selected_D_case, x_fine_mesh, L);
% fD_theor_fine_data = D_func(selected_D_case, x_fine_mesh, L) .* f_func(selected_f_case, x_fine_mesh, L);
% Calculate in bins with first two derivatives
[D_bins, D_prime_bins, D_prime_prime_bins] = D_func(selected_D_case, x_bins_centers, L);
b_bins = sqrt(2 * D_bins);
b_prime_bins = D_prime_bins ./ b_bins;
b_prime_prime_bins = (D_prime_prime_bins .* b_bins - D_prime_bins .* b_prime_bins) ./ b_bins.^2;
% a_bins = f_func(selected_f_case, x_bins_centers, L) / gamma_drag;
% a_theor_fine_data = f_func(selected_f_case, x_fine_mesh, L) / gamma_drag;
b_theor_fine_data = sqrt(2 * D_theor_fine_data);
bb_prime_theor_fine_data = D_grad_theor_fine_data;


fprintf('Processing trajectories...\n');
% Select one bin in the middle to avoid boundary effects
middle_bin = floor(x_bins_number/2);
tic;
parfor trial = 1:trials  % 765
    %% Initialize
    % Initialize the data structure
    data_struct = initialize_data_structure(x_bins_number, fine_mesh_steps_count, conventions_count, lambda_types_count);
% 	data_struct.bl_force = bl_force;
%     data_struct.selected_f_case = selected_f_case;
% 	data_struct.str_force = str_force;
    data_struct.x_bins_number = x_bins_number;
	data_struct.x_bins_centers = x_bins_centers;
    data_struct.x_bins_widths = x_bins_widths;
	data_struct.x_bins_borders = x_bins_borders;
%     data_struct.lambda = trials_lambdas(trial);
    % The following two parameters are the same for the forward and
    % backward calculations because they are just used in the prior
    data_struct.dx_mean_all_bins_all_trials = dx_mean_all_bins_all_trials;
%     data_struct.V = V;
	data_struct.mean_jump_bins_all_trials = mean_jump_bins_all_trials;
    data_struct.x_fine_mesh = x_fine_mesh;
    data_struct.b_theor_fine_data = b_theor_fine_data;
	data_struct.D_theor_fine_data = D_theor_fine_data;
    data_struct.bb_prime_theor_fine_data = bb_prime_theor_fine_data;
%     data_struct.a_theor_fine_data = a_theor_fine_data;
    data_struct.D_theor_data = [D_bins, D_prime_bins, D_prime_prime_bins];
	data_struct.b_theor_data = [b_bins, b_prime_bins, b_prime_prime_bins];
%     data_struct.a_theor_data = a_bins;
%     data_struct.trial_simulation_type = trial_simulation_type;
%     data_struct.trial_first_simulation_type_index = trial_first_simulation_type_index;

    % Ksi-related data
    data_struct.trials_ksi = trials_ksi;
    data_struct.ksi_array = ksi_array;
    data_struct.trials_ksi_type = trials_ksi_type;
    data_struct.trial_first_ksi_type_index = trial_first_ksi_type_index;
    
            
    % Load data for the current trial
    x = trials_x(:, trial)';
    dx = trials_dx(:, trial)';
    
    % Process data with different limits
    data_structs = cell(n_limits_count, 1);
    bin = middle_bin;
    for lim_ind = 1:n_limits_count
        
        n_limit = n_limits(lim_ind);
        cur_data_struct = data_struct;
        
        % Bin data
        [cur_data_struct.n_j, points_in_bins, cur_data_struct.bl_empty_bins] = bin_into_predefined_bins(x, dx, x_bins_borders, n_limit);
        
        % Calculate averages in bins and across bins
        all_kept_dx = [];
        for bin = 1:x_bins_number
            % Individual bins
            cur_data_struct.dx_mean_in_bins(bin) = mean(points_in_bins{bin}(2, :));
            cur_data_struct.V_j(bin) = var(points_in_bins{bin}(2, :));
            cur_data_struct.mean_jump_length_bins(bin) = sqrt(cur_data_struct.V_j(bin));
            
            % Collect data across bins
            all_kept_dx = [all_kept_dx, points_in_bins{bin}(2, :)];
        end
        
        % Average over all bins
        cur_data_struct.dx_mean_all_bins = mean(all_kept_dx);
        cur_data_struct.V_all_bins = var(all_kept_dx);
        
        % Infer MAP b, D and bb'
        cur_data_struct = infer_MAP_b(cur_data_struct);
        
        % Infer force with different conventions
        cur_data_struct = infer_force(cur_data_struct, bin, trial, trials);
        
        % Save data_struct for each n_limit
        data_structs{lim_ind} = cur_data_struct;
    end

    % Save results for this trial
    trials_data{trial} = data_structs;
 
end




data_struct = trials_data{trials};



% Restore the time mesh
t_mesh = (0:N) * t_step;

% Combine predictions from all trials (needed for right parallelization)
trials_MAP_D = zeros(trials, x_bins_number, 4);
trials_MAP_b = zeros(trials, x_bins_number, 4);
trials_MAP_a = zeros(trials, x_bins_number, conventions_count, 4);
trials_MAP_bb_prime_regular_interp = zeros(trials, x_bins_number);
trials_n_j = zeros(trials, x_bins_number);
trials_log_K_L = zeros(trials, x_bins_number, conventions_count);
trials_log_K_G = zeros(trials, conventions_count);
parfor trial = 1:trials
    % D
    trials_MAP_D(trial, :, :) = trials_data{trial}.MAP_D(:, :);
	% b
    trials_MAP_b(trial, :, :) = trials_data{trial}.MAP_b(:, :);
    % fD
    trials_MAP_a(trial, :, :, :) = trials_data{trial}.MAP_a;
    % D grad
    trials_MAP_bb_prime_regular_interp(trial, :) = trials_data{trial}.MAP_bb_prime_regular_interp;
	% n_j
	trials_n_j(trial, :) = trials_data{trial}.n_j;
	% Calculate the Bayes factors K
	fprintf('Calculating the Bayes factor for trial %i/%i\n', trial, trials);
	[log_K_L, log_K_G] = calculate_bayes_factor(trials_data{trial});
	trials_log_K_L(trial, :, :) = log_K_L';
	trials_log_K_G(trial, :) = log_K_G;
end
% Save
data_struct.trials_MAP_D = trials_MAP_D;
data_struct.trials_MAP_b = trials_MAP_b;
data_struct.trials_MAP_a = trials_MAP_a;
data_struct.trials_MAP_bb_prime_regular_interp = trials_MAP_bb_prime_regular_interp;
data_struct.trials_log_K_L = trials_log_K_L;
data_struct.trials_log_K_G = trials_log_K_G;



% % % %% Calculate KS distances for a and b distributions
% % % [trials_a_KS_distance, trials_b_KS_distance] = batch_calculate_KS_distance(trials, x_bins_number, data_struct, trials_data, trials_MAP_a, trials_MAP_b);
% % % 
% % % % Save
% % % data_struct.trials_b_KS_distance = trials_b_KS_distance;



%% Calculate mean for each simulation type separately
%% Also calculate the fail rate for each simulation type
% Initialize arrays
MAP_D_mean = zeros(ksi_count, x_bins_number, 4);
MAP_b_mean = zeros(ksi_count, x_bins_number, 4);
% MAP_a_mean = zeros(lambda_types_count, x_bins_number, conventions_count, 4);
MAP_bb_prime_regular_interp_mean = zeros(ksi_count, x_bins_number);
% UR_b = zeros(lambda_types_count, x_bins_number);
% UR_a = zeros(lambda_types_count, x_bins_number, conventions_count);
% outside_count_a = zeros(lambda_types_count, x_bins_number, conventions_count);
n_j_mean = zeros(ksi_count, x_bins_number);
% b_KS_distance_mean = zeros(lambda_types_count, x_bins_number);
mean_log_K_L = zeros(ksi_count, conventions_count);
std_log_K_L = zeros(ksi_count, conventions_count);
% mean_log_K_G = zeros(ksi_count, conventions_count);

for ksi_ind = 1:ksi_count
    % Mean
    MAP_D_mean(ksi_ind, :, :) = mean(trials_MAP_D(trials_ksi_type == ksi_ind, :, :), 1, 'omitnan' );
	MAP_b_mean(ksi_ind, :, :) = mean(trials_MAP_b(trials_ksi_type == ksi_ind, :, :), 1, 'omitnan' );
%     MAP_a_mean(lambda_type, :, :, :) = mean(trials_MAP_a(trial_simulation_type == lambda_type, :, :, :), 1, 'omitnan' );
    MAP_bb_prime_regular_interp_mean(ksi_ind, :) = mean(trials_MAP_bb_prime_regular_interp(trials_ksi_type == ksi_ind, :), 1, 'omitnan' );
	
	% n_j
	n_j_mean(ksi_ind, :, :) = mean(trials_n_j(trials_ksi_type == ksi_ind, :), 1, 'omitnan' );
        
	% Fail rate b
%     UR_b(lambda_type, :) = mean(double(trials_MAP_b(trial_simulation_type == lambda_type, :, 4) > CONF_LEVEL), 1, 'omitnan' );
    
	% Fail rate a
%     UR_a(lambda_type, :, :) = mean(double(trials_MAP_a(trial_simulation_type == lambda_type, :, :, 4) > CONF_LEVEL), 1, 'omitnan' );
	
% 	% Kolmogorov-Smirnov distance for b
% 	b_KS_distance_mean(lambda_type, :) = mean(trials_b_KS_distance(trial_simulation_type == lambda_type, :), 1, 'omitnan');
	
	% Bayes factors K
% 	mean_log_K_G(lambda_type, :) = mean(trials_log_K_G(trial_simulation_type == lambda_type, :), 1, 'omitnan' );
% 	std_log_K_G(lambda_type, :) = std(trials_log_K_G(trial_simulation_type == lambda_type, :), 1, 'omitnan' );
    for convention = 1:conventions_count
        vals = reshape(trials_log_K_L(trials_ksi_type == ksi_ind, middle_bin, convention),[], 1);
        mean_log_K_L(ksi_ind, convention) = mean(vals, 'omitnan');
        std_log_K_L(ksi_ind, convention) = std(vals, 'omitnan');
    end
end

% Save
data_struct.MAP_D_mean = MAP_D_mean;
data_struct.MAP_b_mean = MAP_b_mean;
% data_struct.MAP_a_mean = MAP_a_mean;
% data_struct.b_KS_distance_mean = b_KS_distance_mean;
% data_struct.b_KS_distance_bin_mean = mean(b_KS_distance_mean, 2);
data_struct.MAP_bb_prime_regular_interp_mean = MAP_bb_prime_regular_interp_mean;
% data_struct.UR_b = UR_b;
% data_struct.UR_a = UR_a;
data_struct.n_j_mean = n_j_mean;
% data_struct.mean_log_K_G = mean_log_K_G;
% data_struct.std_log_K_G = std_log_K_G;
% data_struct.eb_log_K_G = std_log_K_G * sqrt(2) * erfinv(0.95);
data_struct.mean_log_K_L = mean_log_K_L;
data_struct.std_log_K_L = std_log_K_L;
data_struct.eb_log_K_L = std_log_K_L * sqrt(2) * erfinv(0.95);


% %% Calculate mean K_L across x>0 and x<0 excluding the x=0 bin
% zero_bin_number = get_selected_bins_indices(data_struct, 0);
% neg_bins_indices = 1:(zero_bin_number - 1);
% pos_bins_indices = (zero_bin_number + 1):x_bins_number;
% 
% % Average
% mean_log_K_L_avg_neg_x = squeeze(mean(mean_log_K_L(:, neg_bins_indices, :), 2));
% mean_log_K_L_avg_pos_x = squeeze(mean(mean_log_K_L(:, pos_bins_indices, :), 2));
% 
% % Save
% data_struct.mean_log_K_L_avg_neg_x = mean_log_K_L_avg_neg_x;
% data_struct.mean_log_K_L_avg_pos_x = mean_log_K_L_avg_pos_x;




%% Clean up & backup current workspace
trials_x = [];
trials_dx = [];
points_binned = [];
save(strcat('backup_workspace.mat'));
save(strcat(output_data_folder, 'trials_data.mat'), 'data_struct', 'trials_data');

% Print execution time
fprintf('All trajectories processed in %.2f min\n', toc/60);


%% Plot
% plot_article_all(data_struct, trials_data);
fig_count = 1; 
plot_article_local_bayes_factor(data_struct, fig_count, 0);


1;























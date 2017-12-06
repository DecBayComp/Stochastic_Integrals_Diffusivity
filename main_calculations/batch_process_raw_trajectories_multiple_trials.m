set (0, 'DefaultAxesFontSize', 12);


%% Clear the workspace (do not use if want to be able to use bl_reload_trajectories = false)
% clear;

%% Constants
load_constants;

a_ABS_MAX = 10;
D_PRECISION = 1e-5;
b_PRECISION = 1e-3;
D_ABS_MAX = 1;
b_ABS_MAX = 1;
bl_find_marginalized_fD_error_bars = false;		% keep off. A much faster calculation method was implemented
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
if bl_force
% 	input_data_folder = '/home/aserov/Documents/Calculated_data/dilemma_with_force/';
	selected_f_case = enum_force_case;
	str_force = 'with_force';	
else
% 	input_data_folder = '/home/aserov/Documents/Calculated_data/dilemma_no_force/';
	selected_f_case = enum_no_force_case;
	str_force = 'no_force';
end


%% Access trajectories folder and start loading trajectories
% Count the number of csv trajectories in a folder
cur_dir = dir([input_data_folder, '*.csv']);
input_files_count = sum(~[cur_dir.isdir]);
% input_files_count = 40*4;



% Load trajectories
if bl_reload_trajectories
    % Initialize arrays
    trials_x = zeros(N, input_files_count);
    trials_dx = zeros(N, input_files_count);
    trials_lambdas = zeros(1, input_files_count);
	trials_D_case = zeros(1, input_files_count);
	trials_f_case = zeros(1, input_files_count);
    % Load
    parfor file = 1:input_files_count
        filename = sprintf('sim_data_%09i.csv', file);
        fprintf('Loading trajectory from  file %s. Progress: %i/%i\n', filename, file, input_files_count);
        output_full_path = strcat(input_data_folder, filename);

        input_data = dlmread(output_full_path, CSV_DELIMITER);
		
		trials_D_case(file) = input_data(1,1);
		trials_f_case(file) = input_data(1,2);
        trials_lambdas(file) = input_data(2,1);
        trials_x(:, file) = input_data(3:N+2, 1);
        trials_dx(:, file) = input_data(3:N+2,2);
    end;
    input_data = [];
end;



%% Only keep data that corresponds to the currently analyzed condition (force / no force)
trials_x = trials_x(:, trials_f_case == selected_f_case);
trials_dx = trials_dx(:, trials_f_case == selected_f_case);
trials_lambdas = trials_lambdas(trials_f_case == selected_f_case);
trials = sum(trials_f_case == selected_f_case);
trials_data = cell(1, trials);
	


%% Identify and store lambda simulation type
trial_simulation_type = ones(trials, 1) * enum_lambda_rand;
trial_simulation_type (trials_lambdas == 0) = enum_lambda_Ito;
trial_simulation_type (trials_lambdas == 0.5) = enum_lambda_Stratonovich;
trial_simulation_type (trials_lambdas == 1) = enum_lambda_Hanggi;
% Locate first indices of each lambda simulation type
trial_first_simulation_type_index = zeros(lambda_types_count, 1);
for lambda_type = 1:lambda_types_count
    trial_first_simulation_type_index(lambda_type) = find(trial_simulation_type == lambda_type, 1);
end;


% %% Initialize variables (???)
% N_points_total = N * trials;
% x_lambda_all_trials = zeros(lambda_count, N_points_total);
% dx_fwd_lambda_all_trials = zeros(lambda_count, N_points_total);  % Forward jumps
% dx_bck_lambda_all_trials = zeros(lambda_count, N_points_total);  % Backward jumps
% x_lambda = zeros(lambda_count, N);
% dx_fwd_lambda = zeros(lambda_count, N);
% dx_bck_lambda = zeros(lambda_count, N);
% 
% 


%% Identify suitable bin locations based on all points for all trials
[x_bins_borders, x_bins_centers, x_bins_number, x_bins_widths,...
            elements_in_bins_count, variance_in_bins, ~] = select_bins_adaptive_mesh(trials_x(:), trials_dx(:), points_in_bin_avg * trials);

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
fD_theor_fine_data = D_func(selected_D_case, x_fine_mesh, L) .* f_func(selected_f_case, x_fine_mesh, L);
% Calculate in bins with first two derivatives
[D_bins, D_prime_bins, D_prime_prime_bins] = D_func(selected_D_case, x_bins_centers, L);
b_bins = sqrt(2 * D_bins);
b_prime_bins = D_prime_bins ./ b_bins;
b_prime_prime_bins = (D_prime_prime_bins .* b_bins - D_prime_bins .* b_prime_bins) ./ b_bins.^2;
a_bins = f_func(selected_f_case, x_bins_centers, L) / gamma_drag;
a_theor_fine_data = f_func(selected_f_case, x_fine_mesh, L) / gamma_drag;
b_theor_fine_data = sqrt(2 * D_theor_fine_data);
bb_prime_theor_fine_data = D_grad_theor_fine_data;


fprintf('Processing trajectories...\n');
tic;
parfor trial = 1:trials  % 765
    %% Initialize
    % Initialize the data structure
    data_struct = initialize_data_structure(x_bins_number, fine_mesh_steps_count, conventions_count, lambda_types_count);
	data_struct.bl_force = bl_force;
	data_struct.str_force = str_force;
    data_struct.x_bins_number = x_bins_number;
	data_struct.x_bins_centers = x_bins_centers;
    data_struct.x_bins_widths = x_bins_widths;
	data_struct.x_bins_borders = x_bins_borders;
    data_struct.lambda = trials_lambdas(trial);
    % The following two parameters are the same for the forward and
    % backward calculations because they are just used in the prior
    data_struct.dx_mean_all_bins_all_trials = dx_mean_all_bins_all_trials;
%     data_struct.V = V;
	data_struct.mean_jump_bins_all_trials = mean_jump_bins_all_trials;
    data_struct.x_fine_mesh = x_fine_mesh;
    data_struct.b_theor_fine_data = b_theor_fine_data;
	data_struct.D_theor_fine_data = D_theor_fine_data;
    data_struct.bb_prime_theor_fine_data = bb_prime_theor_fine_data;
    data_struct.a_theor_fine_data = a_theor_fine_data;
    data_struct.D_theor_data = [D_bins, D_prime_bins, D_prime_prime_bins];
	data_struct.b_theor_data = [b_bins, b_prime_bins, b_prime_prime_bins];
    data_struct.a_theor_data = a_bins;
    data_struct.trial_simulation_type = trial_simulation_type;
    data_struct.trial_first_simulation_type_index = trial_first_simulation_type_index;
    
    
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
		
		
		%% If for this trial, we have no points in the bin, set flag
		if data_struct.n_j(bin) == 0
			bl_empty_bin = true;
		end
		

        %% Calculate MAP diffusivity D and b
        % Prepare functions
        log_D_function_to_minimze = @(D) bin_D_log_posterior_func (bin, D, t_step, data_struct, 'forward');
		log_b_function_to_minimze = @(b) bin_b_log_posterior_func (bin, b, t_step, data_struct, 'forward');
        % Make an MLE guess
        [mu_n, kappa_n, nu_n, sigma2_n] = get_n_parameters(bin, data_struct, 'forward');
		MLE_guess_D = 2 * nu_n / (nu_n + 2) * sigma2_n / (2 * t_step);
		MLE_guess_b = sqrt(MLE_guess_D);
        % Find confidence intervals if the bin is not empty
		if ~bl_empty_bin
			D_inference = find_confidence_interval(log_D_function_to_minimze, [D_PRECISION, D_ABS_MAX], true, MLE_guess_D, CONF_LEVEL,...
				data_struct.D_theor_data(bin), trial, bin);
			b_inference = find_confidence_interval(log_b_function_to_minimze, [b_PRECISION, b_ABS_MAX], true, MLE_guess_b, CONF_LEVEL,...
				data_struct.b_theor_data(bin), trial, bin);
		else
			D_inference = ones(1, 4) * NaN;
			b_inference = ones(1, 4) * NaN;
		end;
        % Save
        data_struct.MAP_D(bin, :) = D_inference;
		data_struct.MAP_b(bin, :) = b_inference;
		data_struct.bl_empty_bin(bin) = bl_empty_bin;
    end;

    
    %% Regularize bb' gradient based on in non-empty bins
	% To obtain the bb' gradient, provide the function b^2/2 as input
	b_squared_over_2 = data_struct.MAP_b(:, 1).^2 / 2;
% 	if ~bl_empty_bin
	[inferred_MAP_b_squared_over_2_reg, inferred_MAP_bb_prime_reg, inferred_MAP_bb_prime_reg_interpolated, norm_cost, x_grad_mesh] = ...
		regularize_gradient(b_squared_over_2, x_bins_centers, alpha_reg);
% 	else
% 		inferred_MAP_b_squared_over_2_reg = zeros(bins_number, 1) * NaN;
% 		inferred_MAP_bb_prime_reg = zeros(bins_number, 1) * NaN;
% 		inferred_MAP_bb_prime_reg_interpolated = zeros(bins_number, 1) * NaN;
% 		norm_cost = zeros(bins_number, 1) * NaN;
% 		x_grad_mesh = zeros(bins_number, 1) * NaN;
% 	end;
    % Save
    data_struct.MAP_b_regular = sqrt(inferred_MAP_b_squared_over_2_reg * 2);
    data_struct.MAP_bb_prime_regular = inferred_MAP_bb_prime_reg;
    data_struct.MAP_bb_prime_regular_interp = inferred_MAP_bb_prime_reg_interpolated;
    data_struct.x_grad_mesh = x_grad_mesh;

    
    %% Infer forces
    MAP_a = zeros(x_bins_number, conventions_count, 4);
%     MAP_fwd_fD_divine = zeros(4, x_bins_number_saved(MAP_D_grad_regular));
%     MAP_fD_Ito = zeros(4, x_bins_number_saved(MAP_D_grad_regular));
%     MAP_fwd_fD_Stratonovich = zeros(4, x_bins_number_saved(MAP_D_grad_regular));
%     MAP_fD_Hanggi = zeros(4, x_bins_number_saved(MAP_D_grad_regular));
%     MAP_fwd_fD_marginalized = zeros(4, x_bins_number_saved(MAP_D_grad_regular));
    for bin = 1:x_bins_number
        fprintf('Estimating force. Trial: %i/%i. Bin: %i/%i\n', trial, trials, bin, x_bins_number);
        % Initialize
        [mu_n, kappa_n, nu_n, sigma2_n] = get_n_parameters(bin, data_struct, 'forward');
        bb_prime = inferred_MAP_bb_prime_reg_interpolated(bin);
		bl_empty_bin = data_struct.bl_empty_bin(bin);
		
		
		
		%% Skip if bin is empty
		if bl_empty_bin
			% Set to NaN
			a_divine_inference = ones(1, 4) * NaN;
			a_Ito_inference = ones(1, 4) * NaN;
			a_Stratonovich_inference = ones(1, 4) * NaN;
			a_Hanggi_inference = ones(1, 4) * NaN;
			a_MLE_marginalized = ones(1, 4) * NaN;
			
			% Save
			MAP_a(bin, enum_conv_divine, :) = a_divine_inference;
			MAP_a(bin, enum_conv_Ito, :) = a_Ito_inference;
			MAP_a(bin, enum_conv_Stratonovich, :) = a_Stratonovich_inference;
			MAP_a(bin, enum_conv_Hanggi, :) = a_Hanggi_inference;
			MAP_a(bin, enum_conv_marginalized, :) = a_MLE_marginalized;
			
			display('This bin is empty. Skipping');
			continue;
		end;
		


        %% Divine force estimate
        % Prepare function
        log_function_to_minimze = @(a) bin_a_divine_inference_log_posterior_func(data_struct, trials_lambdas(trial), bin, a, bb_prime, 'forward');
        % Make an MLE guess
        lambda = data_struct.lambda;
        MLE_guess = mu_n / t_step - lambda * bb_prime;
        % Find confidence intervals if bin not empty
% 		if ~
        a_divine_inference = find_confidence_interval(log_function_to_minimze, [- a_ABS_MAX, a_ABS_MAX], true, MLE_guess,...
            CONF_LEVEL, data_struct.a_theor_data(bin), trial, bin);
        % Save
        MAP_a(bin, enum_conv_divine, :) = a_divine_inference;

		
		
        %% Simple Ito force estimate
        % Prepare function
        function_to_minimze = @(a) bin_a_log_posterior_func (data_struct, bin, a, 'forward');
        % Make an MLE guess
        MLE_guess = mu_n / t_step;
        % Find confidence intervals
        a_Ito_inference = find_confidence_interval(function_to_minimze, [- a_ABS_MAX, a_ABS_MAX], true, MLE_guess,...
            CONF_LEVEL, data_struct.a_theor_data(bin), trial, bin);
        % Save
        MAP_a(bin, enum_conv_Ito, :) = a_Ito_inference;
		
		

        %% Simple Stratonovich force estimate (through Ito)
        % Prepare function
        function_to_minimze = @(a) bin_a_simple_Stratonovich_log_posterior_func(data_struct, bin, a, bb_prime, 'forward');
        % Make an MLE guess
        lambda = 1/2;
        MLE_guess = mu_n / t_step - lambda * bb_prime;
        % Find confidence intervals
        a_Stratonovich_inference = find_confidence_interval(function_to_minimze, [- a_ABS_MAX, a_ABS_MAX], true, MLE_guess,...
            CONF_LEVEL, data_struct.a_theor_data(bin), trial, bin);
        % Save
        MAP_a(bin, enum_conv_Stratonovich, :) = a_Stratonovich_inference;
		
		

        %% Simple Hanggi force estimate (through Ito)
        % Prepare function
        function_to_minimze = @(a) bin_a_simple_Hanggi_log_posterior_func(data_struct, bin, a, bb_prime, 'forward');
        % Make an MLE guess
        lambda = 1;
        MLE_guess = mu_n / t_step - lambda * bb_prime;
        % Find confidence intervals
        a_Hanggi_inference = find_confidence_interval(function_to_minimze, [- a_ABS_MAX, a_ABS_MAX], true, MLE_guess,...
            CONF_LEVEL, data_struct.a_theor_data(bin), trial, bin);
        % Save
        MAP_a(bin, enum_conv_Hanggi, :) = a_Hanggi_inference;
		
		

        %% Marginalized force estimate (through Ito)
        % Prepare function
        function_to_minimze = @(a) bin_a_lambda_marginalized_log_posterior_func(data_struct, bin, a, bb_prime, 'forward');
        % Make a guess
        lambda = 0.5;
        MLE_guess = mu_n / t_step - lambda * bb_prime;
        % Find confidence intervals
        a_MLE_marginalized = find_confidence_interval(function_to_minimze, [- a_ABS_MAX, a_ABS_MAX], bl_find_marginalized_fD_error_bars,...
            MLE_guess, CONF_LEVEL, data_struct.a_theor_data(bin), trial, bin);
        % Save
        MAP_a(bin, enum_conv_marginalized, :) = a_MLE_marginalized;

    end;
    
    % Save all to the data structure 
    data_struct.MAP_a = MAP_a;
    
    %% Save results for this trial
    trials_data{trial} = data_struct;
 
end;
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
end;
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
MAP_D_mean = zeros(lambda_types_count, x_bins_number, 4);
MAP_b_mean = zeros(lambda_types_count, x_bins_number, 4);
MAP_a_mean = zeros(lambda_types_count, x_bins_number, conventions_count, 4);
MAP_bb_prime_regular_interp_mean = zeros(lambda_types_count, x_bins_number);
UR_b = zeros(lambda_types_count, x_bins_number);
UR_a = zeros(lambda_types_count, x_bins_number, conventions_count);
outside_count_a = zeros(lambda_types_count, x_bins_number, conventions_count);
n_j_mean = zeros(lambda_types_count, x_bins_number);
b_KS_distance_mean = zeros(lambda_types_count, x_bins_number);
mean_log_K_L = zeros(lambda_types_count, x_bins_number, conventions_count);
mean_log_K_G = zeros(lambda_types_count, conventions_count);

for lambda_type = 1:lambda_types_count
    % Mean
    MAP_D_mean(lambda_type, :, :) = mean(trials_MAP_D(trial_simulation_type == lambda_type, :, :), 1, 'omitnan' );
	MAP_b_mean(lambda_type, :, :) = mean(trials_MAP_b(trial_simulation_type == lambda_type, :, :), 1, 'omitnan' );
    MAP_a_mean(lambda_type, :, :, :) = mean(trials_MAP_a(trial_simulation_type == lambda_type, :, :, :), 1, 'omitnan' );
    MAP_bb_prime_regular_interp_mean(lambda_type, :) = mean(trials_MAP_bb_prime_regular_interp(trial_simulation_type == lambda_type, :), 1, 'omitnan' );
	
	% n_j
	n_j_mean(lambda_type, :, :) = mean(trials_n_j(trial_simulation_type == lambda_type, :), 1, 'omitnan' );
        
	% Fail rate b
    UR_b(lambda_type, :) = mean(double(trials_MAP_b(trial_simulation_type == lambda_type, :, 4) > CONF_LEVEL), 1, 'omitnan' );
    
	% Fail rate a
    UR_a(lambda_type, :, :) = mean(double(trials_MAP_a(trial_simulation_type == lambda_type, :, :, 4) > CONF_LEVEL), 1, 'omitnan' );
	
% 	% Kolmogorov-Smirnov distance for b
% 	b_KS_distance_mean(lambda_type, :) = mean(trials_b_KS_distance(trial_simulation_type == lambda_type, :), 1, 'omitnan');
	
	% Bayes factors K
	mean_log_K_L(lambda_type, :, :) = mean(trials_log_K_L(trial_simulation_type == lambda_type, :, :), 1, 'omitnan');
	mean_log_K_G(lambda_type, :) = mean(trials_log_K_G(trial_simulation_type == lambda_type, :), 1, 'omitnan' );
end;

% Save
data_struct.MAP_D_mean = MAP_D_mean;
data_struct.MAP_b_mean = MAP_b_mean;
data_struct.MAP_a_mean = MAP_a_mean;
% data_struct.b_KS_distance_mean = b_KS_distance_mean;
% data_struct.b_KS_distance_bin_mean = mean(b_KS_distance_mean, 2);
data_struct.MAP_bb_prime_regular_interp_mean = MAP_bb_prime_regular_interp_mean;
data_struct.UR_b = UR_b;
data_struct.UR_a = UR_a;
data_struct.n_j_mean = n_j_mean;
data_struct.mean_log_K_L = mean_log_K_L;
data_struct.mean_log_K_G = mean_log_K_G;



%% Calculate mean and average fail rate of each inference convention
% The mean is performed over bins and is weighted with bins size

% Identify the best explored period
x_left = (1/2 + 2*1)/w;
x_right = (1/2 + 2*2)/w;

% Filter indices from one best period
indices = data_struct.x_bins_centers >= x_left & data_struct.x_bins_centers <= x_right;
bin_widths = data_struct.x_bins_centers(indices);
norm_bin_widths = bin_widths / sum(bin_widths);

% UR_fD indices: (lambda_types_count, x_bins_number, conventions_count)
UR_b_bin_mean = zeros(lambda_types_count, 1);
UR_a_bin_mean = zeros(lambda_types_count, conventions_count);
A = permute(data_struct.UR_a, [1, 3, 2]);
for lambda_type = 1:lambda_types_count
    UR_b_bin_mean(lambda_type) = data_struct.UR_b(lambda_type, indices) * norm_bin_widths;
    UR_a_bin_mean(lambda_type, :) = squeeze(A(lambda_type, :, indices)) * norm_bin_widths;
end;
UR_b_bin_max = max(data_struct.UR_b(:, indices), [], 2);
UR_a_bin_max = squeeze(max(data_struct.UR_a(:, indices, :), [], 2));

% Save
data_struct.UR_b_bin_mean = UR_b_bin_mean;
data_struct.UR_b_bin_max = UR_b_bin_max;
data_struct.UR_a_bin_mean = UR_a_bin_mean;
data_struct.UR_a_bin_max = UR_a_bin_max;




%% Clean up & backup current workspace
trials_x = [];
trials_dx = [];
points_binned = [];
save(strcat('backup_workspace_', str_force, '.mat'));
save(strcat(output_data_folder, 'trials_data.mat'), 'data_struct', 'trials_data');

% Print execution time
fprintf('All trajectories processed in %.2f min\n', toc/60);


%% Plot
plot_article_all(data_struct, trials_data);


1;























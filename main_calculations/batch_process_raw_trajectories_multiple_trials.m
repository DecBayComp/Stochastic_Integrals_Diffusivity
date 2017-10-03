set (0, 'DefaultAxesFontSize', 12);


%% Clear the workspace
clear;

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


%% Access trajectories folder and start loading trajectories
% Count the number of csv trajectories in a folder
cur_dir = dir([input_data_folder, '*.csv']);
trials = sum(~[cur_dir.isdir]);
% trials = 11 * 10;


% Load trajectories
if bl_reload_trajectories
    % Initialize arrays
    trials_data = cell(1, trials);
    trials_x = zeros(N, trials);
    trials_dx = zeros(N, trials);
    trials_lambdas = zeros(1, trials);
    % Load
    parfor trial = 1:trials
        filename = sprintf('sim_data_%09i.csv', trial);
        fprintf('Loading trajectory from  file %s. Progress: %i/%i\n', filename, trial, trials);
        output_full_path = strcat(input_data_folder, filename);

        input_data = dlmread(output_full_path, CSV_DELIMITER);

        trials_lambdas(trial) = input_data(1);
        trials_x(:, trial) = input_data(2:N+1);
        trials_dx(:, trial) = input_data(3:N+2) - input_data(2:N+1);
    end;
    input_data = [];
end;

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
            elements_in_bins_count, ~] = select_bins_adaptive_mesh(trials_x(:), trials_dx(:), points_in_bin_avg * trials);
% Estimate dx_Mean and V used for prior only. Average over everything
dx_Mean = mean(trials_dx(:));
V = mean(trials_dx(:).^2) - dx_Mean^2;


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
a_bins = f_func(selected_f_case, x_bins_centers, L) / gamma_drag;
a_theor_fine_data = f_func(selected_f_case, x_fine_mesh, L) / gamma_drag;
b_theor_fine_data = sqrt(2 * D_theor_fine_data);
bb_prime_theor_fine_data = D_grad_theor_fine_data;


fprintf('Processing trajectories...\n');
tic;
parfor trial = 1:trials
    %% Initialize
    % Initialize the data structure
    data_struct = initialize_data_structure(x_bins_number, fine_mesh_steps_count, conventions_count);
    data_struct.x_bins_centers = x_bins_centers;
    data_struct.x_bins_widths = x_bins_widths;
    data_struct.lambda = trials_lambdas(trial);
    % The following two parameters are the same for the forward and
    % backward calculations because they are just used in the prior
    data_struct.dx_Mean = dx_Mean;
    data_struct.V = V;
    data_struct.x_fine_mesh = x_fine_mesh;
    data_struct.b_theor_fine_data = b_theor_fine_data;
    data_struct.bb_prime_theor_fine_data = bb_prime_theor_fine_data;
    data_struct.a_theor_fine_data = a_theor_fine_data;
    data_struct.D_theor_data = [D_bins; D_prime_bins; D_prime_prime_bins];
	data_struct.b_theor_data(:, 1) = b_bins;
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
        elements_in_bins_count(bin) = sum((sorted_data(1, :) > x_bins_borders(bin, 1)) & (sorted_data(1, :) <= x_bins_borders(bin, 2)));
        counter_end = counter_start + elements_in_bins_count(bin) - 1;
        points_in_bins{bin} = sorted_data(:, counter_start : counter_end);
        counter_start = counter_end + 1;

        %% Save calculated values to the data structure
        data_struct.dx_mean_in_bins(bin) = mean(points_in_bins{bin}(2, :));
%         data_struct.dx_bck_mean_in_bins_saved{l_ind}(bin) = mean(points_in_bins{bin}(3, :));
        data_struct.n_j(bin) = elements_in_bins_count(bin);
        data_struct.V_j(bin) = mean(points_in_bins{bin}(2, :).^2) - mean(points_in_bins{bin}(2, :))^2;
%         data_struct.V_bck_j{l_ind}(bin) = mean(points_in_bins{bin}(3, :).^2) - mean(points_in_bins{bin}(3, :))^2;

        %% Calculate MAP diffusivity
        [mu_n, kappa_n, nu_n, sigma2_n] = get_n_parameters(bin, data_struct, 'forward');
        % Prepare function
        log_function_to_minimze = @(b) bin_b_log_posterior_func (bin, b, t_step, data_struct, 'forward');
        % Make an MLE guess
        MLE_guess = sqrt(2 * nu_n / (nu_n + 2) * sigma2_n / (2 * t_step));
        % Find confidence intervals
        b_inference = find_confidence_interval(log_function_to_minimze, [b_PRECISION, b_ABS_MAX], true, MLE_guess, CONF_LEVEL,...
            data_struct.b_theor_data(bin), trial, bin);
        % Save
        data_struct.MAP_b(bin, :) = b_inference;
    end;

    
    %% Regularize bb' gradient based on forward calculations
	% To obtain the bb' gradient, provide the function b^2/2 as input
	b_squared_over_2 = data_struct.MAP_b(:, 1).^2 / 2;
    [inferred_MAP_b_squared_over_2_reg, inferred_MAP_bb_prime_reg, inferred_MAP_bb_prime_reg_interpolated, norm_cost, x_grad_mesh] = ...
        regularize_gradient(b_squared_over_2, x_bins_centers, alpha_reg);
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


        %% Divine force estimate
        % Prepare function
        log_function_to_minimze = @(a) bin_a_divine_inference_log_posterior_func(data_struct, trials_lambdas(trial), bin, a, bb_prime, 'forward');
        % Make an MLE guess
        lambda = data_struct.lambda;
        MLE_guess = mu_n / t_step - lambda * bb_prime;
        % Find confidence intervals
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
data_struct = trials_data{end};

% Restore the time mesh
t_mesh = (0:N) * t_step;

% Combine predictions from all trials (needed for right parallelization)
trials_MAP_b = zeros(trials, x_bins_number, 4);
trials_MAP_a = zeros(trials, x_bins_number, conventions_count, 4);
trials_MAP_bb_prime_regular_interp = zeros(trials, x_bins_number);
for trial = 1:trials
    % b
    trials_MAP_b(trial, :, :) = trials_data{trial}.MAP_b(:, :);
    % fD
    trials_MAP_a(trial, :, :, :) = trials_data{trial}.MAP_a;
    % D grad
    trials_MAP_bb_prime_regular_interp(trial, :) = trials_data{trial}.MAP_bb_prime_regular_interp;
end;
% Save
data_struct.trials_MAP_b = trials_MAP_b;
data_struct.trials_MAP_a = trials_MAP_a;
data_struct.trials_MAP_bb_prime_regular_interp = trials_MAP_bb_prime_regular_interp;



%% Calculate mean for each simulation type separately
%% Also calculate the fail rate for each simulation type
MAP_b_mean = zeros(lambda_types_count, x_bins_number, 4);
MAP_a_mean = zeros(lambda_types_count, x_bins_number, conventions_count, 4);
MAP_bb_prime_regular_interp_mean = zeros(lambda_types_count, x_bins_number);
UR_b = zeros(lambda_types_count, x_bins_number);
UR_a = zeros(lambda_types_count, x_bins_number, conventions_count);
outside_count_a = zeros(lambda_types_count, x_bins_number, conventions_count);
for lambda_type = 1:lambda_types_count
    % Mean
    MAP_b_mean(lambda_type, :, :) = mean(trials_MAP_b(trial_simulation_type == lambda_type, :, :), 1);
    MAP_a_mean(lambda_type, :, :, :) = mean(trials_MAP_a(trial_simulation_type == lambda_type, :, :, :), 1);
    MAP_bb_prime_regular_interp_mean(lambda_type, :) = mean(trials_MAP_bb_prime_regular_interp(trial_simulation_type == lambda_type, :), 1);
    % Fail rate
    % b
    UR_b(lambda_type, :) = mean(double(trials_MAP_b(trial_simulation_type == lambda_type, :, 4) > CONF_LEVEL), 1);
    % a
    UR_a(lambda_type, :, :) = mean(double(trials_MAP_a(trial_simulation_type == lambda_type, :, :, 4) > CONF_LEVEL), 1);
end;
% Save
data_struct.MAP_b_mean = MAP_b_mean;
data_struct.MAP_a_mean = MAP_a_mean;
data_struct.MAP_bb_prime_regular_interp_mean = MAP_bb_prime_regular_interp_mean;
data_struct.UR_b = UR_b;
data_struct.UR_a = UR_a;


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


% Print execution time
fprintf('All trajectories processed in %.2f min\n', toc/60);


%% Clean up & backup current workspace
trials_x = [];
trials_dx = [];
save('backup_workspace.mat');
save(strcat(output_data_folder, 'trials_data.mat'), 'data_struct', 'trials_data');


%% Plot
plot_article_all(data_struct, trials_data);


1;























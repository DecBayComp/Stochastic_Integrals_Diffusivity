set (0, 'DefaultAxesFontSize', 12);


%% Globals
global dx_Mean;
global V;


%% Constants
load_constants;
trials = 0;
fD_ABS_MAX = 4 * 10;
D_PRECISION = 1e-5;
D_ABS_MAX = 1;
bl_find_marginalized_fD_error_bars = false;


%% Initialize
low_conf_percentile = (1 - CONF_LEVEL)/2;
high_conf_percentile = 1 - low_conf_percentile;
fD_min = [];
fD_max = [];
fD_step = [];
fD_mesh = [];
fD_pdf_plot_data = [];
pdf_norm = [];

% Note that N and N_points_total is the number of points and not the number
% of jumps!
N_points_total = N * trials;
trials_data = cell(1, trials);

%% Access trajectories folder and start loading trajectories
% Set first filename
filename = sprintf('sim_data', selected_D_case, selected_f_case, lambda_array(l_ind));

% fprintf('Processing case. D: %i/%i. f: %i/%i. Lambda: %i/%i\n',...


% for D_case_number = 1:max_D_case_number
%     for f_case_number = 1:max_f_case_number
% selected_D_case = selected_D_case;
% selected_f_case = selected_f_case;
%% Initialize variables
x_lambda_all_trials = zeros(lambda_count, N_points_total);
dx_fwd_lambda_all_trials = zeros(lambda_count, N_points_total);  % Forward jumps
dx_bck_lambda_all_trials = zeros(lambda_count, N_points_total);  % Backward jumps
x_lambda = zeros(lambda_count, N);
dx_fwd_lambda = zeros(lambda_count, N);
dx_bck_lambda = zeros(lambda_count, N);
% Load the data for all trials for all lambda values
parfor l_ind = 1:lambda_count
    fprintf('Processing case. D: %i/%i. f: %i/%i. Lambda: %i/%i\n',...
        selected_D_case, max_D_case_number, selected_f_case, max_f_case_number, l_ind, lambda_count);
    %% Load CSV file with x and dx data
    filename = sprintf('D_%i_f_%i_lambda_%.2f_trajectory.csv', selected_D_case, selected_f_case, lambda_array(l_ind));
    output_full_path = strcat(input_data_folder, filename);

    input_data = dlmread(output_full_path, CSV_DELIMITER);
%             if N_total == 0
%                 N_total = size(input_data, 1);
%             end;
    % Separate data into x and dx for both forward and backward jumps
    x_lambda_all_trials(l_ind,:) = input_data(1:N_points_total, 1);
    dx_fwd_lambda_all_trials(l_ind, :) = [input_data(2:N_points_total, 1) - input_data(1:N_points_total-1, 1); 0];
    dx_bck_lambda_all_trials(l_ind, :) = [0; input_data(2:N_points_total, 1) - input_data(1:N_points_total-1, 1)];

end;
input_data = [];


1;


%% Identify suitable bin location for each lambda for the first trial 
%% Using adaptive mesh and minimum number of points per bin
% Copy data
x_lambda(:, 1:N) = x_lambda_all_trials(:, 1:N);
dx_fwd_lambda(:, 1:N) = dx_fwd_lambda_all_trials(:, 1:N);
% Bin
bin_and_infer_bin_distributions;
% This calculates:
% x_bins_borders_saved, x_bins_centers_saved, x_bins_number_saved, x_bins_widths_saved
%
% These bins are the same for both forward and backward jumps (Ito and
% Hanggi calculations) to allow comparison between direct Ito and Hanggi
% inference



%% Once the bins are idetnified, bin each trial with the same frozen bin borders
% Initialize
MAP_D_mean = cell(1, lambda_count);
MAP_fwd_fD_Stratonovich_mean = cell(1, lambda_count);
MAP_fD_Ito_mean = cell(1, lambda_count);
MAP_fwd_fD_divine_mean = cell(1, lambda_count);
MAP_fwd_fD_marginalized_mean = cell(1, lambda_count);
MAP_fD_Hanggi_mean = cell(1, lambda_count);
s1 = 4;
for l_ind = 1:lambda_count
    MAP_D_mean{l_ind} = zeros(s1, x_bins_number_saved(l_ind));
    MAP_fwd_fD_Stratonovich_mean{l_ind} = zeros(s1, x_bins_number_saved(l_ind));
    MAP_fD_Ito_mean{l_ind} = zeros(s1, x_bins_number_saved(l_ind));
    MAP_fwd_fD_divine_mean{l_ind} = zeros(s1, x_bins_number_saved(l_ind));
    MAP_fwd_fD_marginalized_mean{l_ind} = zeros(s1, x_bins_number_saved(l_ind));
    MAP_fD_Hanggi_mean{l_ind} = zeros(s1, x_bins_number_saved(l_ind));
end;


for trial = 1:trials

    %% Initialize
    inferred_MAP_D = cell(1, lambda_count);
    for l_ind = 1:lambda_count
        inferred_MAP_D{l_ind} = zeros(1, x_bins_number_saved(l_ind));
    end;
    % Initialize the data structure
    data_struct = initialize_data_structure(lambda_count, x_bins_number_saved, fine_mesh_steps_count);
    % V and dx_Mean are mean parameters used for prior construction. They will
    % be the same for all trials
    data_struct.x_bins_centers_saved = x_bins_centers_saved;
    data_struct.x_bins_widths_saved = x_bins_widths_saved;
    % The following two parameters are the same for the forward and
    % backward calculations because they are just used in the prior
    data_struct.dx_Mean = dx_Mean;
    data_struct.V = V;
    % Prepare the fine mesh
    x_theor_min = x_min;
    x_theor_max = x_max;
    x_theor_step = (x_theor_max - x_theor_min) / (fine_mesh_steps_count - 1);
    data_struct.x_fine_mesh = x_theor_min:x_theor_step:x_theor_max;
    %% Calculate true D and fD values on a fine mesh and in bins
    % Calculate on the fine mesh
    [data_struct.D_theor_fine_data, data_struct.D_grad_theor_fine_data] = D_func(selected_D_case, data_struct.x_fine_mesh, L);
    data_struct.fD_theor_fine_data = D_func(selected_D_case, data_struct.x_fine_mesh, L) .* f_func(selected_f_case, data_struct.x_fine_mesh, L);
    % Calculate in bins
    for l_ind = 1:lambda_count
        % Getting D and its first two derivatives
        [r1, r2, r3] = D_func(selected_D_case, data_struct.x_bins_centers_saved{l_ind}, L);
        data_struct.D_theor_data{l_ind} = [r1'; r2'; r3'];
        data_struct.fD_theor_data{l_ind} = D_func(selected_D_case, data_struct.x_bins_centers_saved{l_ind}, L)...
            .* f_func(selected_f_case, data_struct.x_bins_centers_saved{l_ind}, L);
    end;
    
    
    
    for l_ind = 1:lambda_count
        tic;
        %% Bin the points for the current trial
        x_data(1, 1:N) = x_lambda_all_trials(l_ind, (1:N) + (trial - 1) * N);
        dx_fwd_data(1, 1:N) = dx_fwd_lambda_all_trials(l_ind, (1:N) + (trial - 1) * N);
        dx_bck_data(1, 1:N) = dx_bck_lambda_all_trials(l_ind, (1:N) + (trial - 1) * N);
        % First sort data
        [~, sorted_indices] = sort(x_data, 'ascend');
        sorted_data = [x_data(sorted_indices); dx_fwd_data(sorted_indices); dx_bck_data(sorted_indices)];
        % Initialize binned data cell array
        points_in_bins = cell(1, x_bins_number_saved(l_ind));
        elements_in_bins_count = zeros(1, x_bins_number_saved(l_ind));
        x_bins_borders = x_bins_borders_saved{l_ind};
        x_bins_centers = x_bins_centers_saved{l_ind};
        % Bin data
        counter_start = 1;
        for bin = 1:x_bins_number_saved(l_ind)
            elements_in_bins_count(bin) = sum((sorted_data(1, :) > x_bins_borders(bin, 1)) & (sorted_data(1, :) <= x_bins_borders(bin, 2)));
            counter_end = counter_start + elements_in_bins_count(bin) - 1;
            points_in_bins{bin} = sorted_data(:, counter_start : counter_end);
            counter_start = counter_end + 1;

            %% Save to the data structure
            data_struct.dx_fwd_mean_in_bins_saved{l_ind}(bin) = mean(points_in_bins{bin}(2, :));
            data_struct.dx_bck_mean_in_bins_saved{l_ind}(bin) = mean(points_in_bins{bin}(3, :));
            data_struct.n_j{l_ind}(bin) = elements_in_bins_count(bin);
            data_struct.V_fwd_j{l_ind}(bin) = mean(points_in_bins{bin}(2, :).^2) - mean(points_in_bins{bin}(2, :))^2;
            data_struct.V_bck_j{l_ind}(bin) = mean(points_in_bins{bin}(3, :).^2) - mean(points_in_bins{bin}(3, :))^2;

            %% Calculate MAP diffusivity
            [mu_fwd_n, kappa_n, nu_n, sigma2_n] = get_n_parameters(l_ind, bin, data_struct, 'forward');
             % Prepare function
            log_function_to_minimze = @(D) bin_D_log_posterior_func (l_ind, bin, D, t_step, data_struct, 'forward');
            % Make an MLE guess
            MLE_guess = nu_n / (nu_n + 2) * sigma2_n / (2 * t_step);
            % Find confidence intervals
            D_inference = find_confidence_interval(log_function_to_minimze, [D_PRECISION, D_ABS_MAX], true, MLE_guess, CONF_LEVEL,...
                data_struct.D_theor_data{l_ind}(bin));
            % Save
            data_struct.MAP_fwd_D{l_ind}(:, bin) = D_inference;
        end;

        % Regularize gradient based on forward calculations
        [inferred_MAP_fwd_D_reg, inferred_MAP_fwd_D_grad_reg, inferred_MAP_fwd_D_grad_reg_interpolated, norm_cost, x_grad_mesh] = ...
            regularize_gradient(data_struct.MAP_fwd_D{l_ind}(1, :), x_bins_centers, alpha_reg);
        % Save
        data_struct.MAP_fwd_D_regular{l_ind} = inferred_MAP_fwd_D_reg;
        data_struct.MAP_fwd_D_grad_regular{l_ind} = inferred_MAP_fwd_D_grad_reg;
        data_struct.MAP_fwd_D_grad_regular_interp{l_ind} = inferred_MAP_fwd_D_grad_reg_interpolated;
        data_struct.x_grad_mesh{l_ind} = x_grad_mesh;
        



        %% Infer forces
        MAP_fwd_fD_divine = zeros(4, x_bins_number_saved(l_ind));
        MAP_fD_Ito = zeros(4, x_bins_number_saved(l_ind));
        MAP_fwd_fD_Stratonovich = zeros(4, x_bins_number_saved(l_ind));
        MAP_fD_Hanggi = zeros(4, x_bins_number_saved(l_ind));
        MAP_fwd_fD_marginalized = zeros(4, x_bins_number_saved(l_ind));
        parfor bin = 1:x_bins_number_saved(l_ind)
            fprintf('Estimating force. Trial: %i/%i. Lambda: %i/%i. Bin: %i/%i\n',...
                trial, trials, l_ind, lambda_count, bin, x_bins_number);
            % Initialize
            [mu_fwd_n, kappa_n, nu_n, sigma2_n] = get_n_parameters(l_ind, bin, data_struct, 'forward');
            D_fwd_grad = inferred_MAP_fwd_D_grad_reg_interpolated(bin);
            
           
            %% Divine force estimate
            % Prepare function
            log_function_to_minimze = @(fD) bin_fD_divine_inference_log_posterior_func(data_struct, ...
                l_ind, bin, fD, inferred_MAP_fwd_D_grad_reg_interpolated(bin), 'forward');
            % Make an MLE guess
            MLE_guess = (mu_fwd_n / t_step - (lambda_array(l_ind)) * D_fwd_grad) * kBT;
            % Find confidence intervals
            fD_divine_inference = find_confidence_interval(log_function_to_minimze, [- fD_ABS_MAX, fD_ABS_MAX], true, MLE_guess,...
                CONF_LEVEL, data_struct.fD_theor_data{l_ind}(bin));
            % Save
            MAP_fwd_fD_divine(:, bin) = fD_divine_inference;
            
            %% Simple Ito force estimate
            % Prepare function
            function_to_minimze = @(fD) bin_fD_log_posterior_func (data_struct, l_ind, bin, fD, 'forward');
            % Make an MLE guess
            MLE_guess = (mu_fwd_n / t_step) * kBT;
            % Find confidence intervals
            fD_Ito_inference = find_confidence_interval(function_to_minimze, [- fD_ABS_MAX, fD_ABS_MAX], true, MLE_guess,...
                CONF_LEVEL, data_struct.fD_theor_data{l_ind}(bin));
    %         fprintf('FINISHED: Calculating ITO inference\n');
            % Save
            MAP_fD_Ito(:, bin) = fD_Ito_inference;

            %% Simple Stratonovich force estimate (through Ito)
            % Prepare function
            function_to_minimze = @(fD) bin_fD_simple_Stratonovich_log_posterior_func(data_struct, ...
                l_ind, bin, fD, inferred_MAP_fwd_D_grad_reg_interpolated(bin), 'forward');
            % Make an MLE guess
            l = 1/2;
            MLE_guess = (mu_fwd_n / t_step - l * D_fwd_grad) * kBT;
            % Find confidence intervals
            fD_Stratonovich_inference = find_confidence_interval(function_to_minimze, [- fD_ABS_MAX, fD_ABS_MAX], true, MLE_guess,...
                CONF_LEVEL, data_struct.fD_theor_data{l_ind}(bin));
            % Save
            MAP_fwd_fD_Stratonovich(:, bin) = fD_Stratonovich_inference;
            
            %% Simple Hanggi force estimate (through Ito)
            % Prepare function
            function_to_minimze = @(fD) bin_fD_simple_Hanggi_log_posterior_func(data_struct, ...
                l_ind, bin, fD, inferred_MAP_fwd_D_grad_reg_interpolated(bin), 'forward');
            % Make an MLE guess
            l = 1;
            MLE_guess = (mu_fwd_n / t_step - l * D_fwd_grad) * kBT;
            % Find confidence intervals
            fD_Hanggi_inference = find_confidence_interval(function_to_minimze, [- fD_ABS_MAX, fD_ABS_MAX], true, MLE_guess,...
                CONF_LEVEL, data_struct.fD_theor_data{l_ind}(bin));
            % Save
            MAP_fD_Hanggi(:, bin) = fD_Hanggi_inference;

            %% Marginalized force estimate (through Ito)
            % Prepare function
            function_to_minimze = @(fD) bin_fD_lambda_marginalized_log_posterior_func(data_struct, ...
                l_ind, bin, fD, inferred_MAP_fwd_D_grad_reg_interpolated(bin), 'forward');
            % Make a guess
            [mu_fwd_n, ~, ~, ~] = get_n_parameters(l_ind, bin, data_struct, 'forward');
            D_fwd_grad = inferred_MAP_fwd_D_grad_reg_interpolated(bin);
            MLE_guess = (mu_fwd_n / t_step - (1/2) * D_fwd_grad) * kBT;
            % Find confidence intervals
            fD_MLE_marginalized = find_confidence_interval(function_to_minimze, [- fD_ABS_MAX, fD_ABS_MAX], bl_find_marginalized_fD_error_bars,...
                MLE_guess, CONF_LEVEL, data_struct.fD_theor_data{l_ind}(bin));
            % Save
            MAP_fwd_fD_marginalized(:, bin) = fD_MLE_marginalized;
            
%             %% Direct Hänggi force estimate
%             % Prepare function
%             function_to_minimze = @(fD) bin_fD_log_posterior_func (data_struct, l_ind, bin, fD, 'backward');
%             % Make an MLE guess
%             [mu_bck_n, ~, ~, ~] = get_n_parameters(l_ind, bin, data_struct, 'backward');
%             MLE_guess = (mu_bck_n / t_step) * kBT;
%             % Find confidence intervals
%             fD_bck_Hanggi_inference = find_confidence_interval(function_to_minimze, [- fD_ABS_MAX, fD_ABS_MAX], true, MLE_guess);
%     %         fprintf('FINISHED: Calculating ITO inference\n');
%             % Save
%             data_struct.MAP_bck_fD_Hanggi{l_ind}(:, bin) = fD_bck_Hanggi_inference;


        end;
        % Save to the data structure
        data_struct.MAP_fwd_fD_divine{l_ind} = MAP_fwd_fD_divine;
        data_struct.MAP_fD_Ito{l_ind} = MAP_fD_Ito;
        data_struct.MAP_fwd_fD_Stratonovich{l_ind} = MAP_fwd_fD_Stratonovich;
        data_struct.MAP_fD_Hanggi{l_ind} = MAP_fD_Hanggi;
        data_struct.MAP_fwd_fD_marginalized{l_ind} = MAP_fwd_fD_marginalized;
        
        % Calculate mean
        MAP_D_mean{l_ind} = MAP_D_mean{l_ind} + data_struct.MAP_fwd_D{l_ind};
%         MAP_bck_fD_Hanggi_mean{l_ind} = MAP_bck_fD_Hanggi_mean{l_ind} + data_struct.MAP_bck_fD_Hanggi{l_ind};
        MAP_fwd_fD_Stratonovich_mean{l_ind} = MAP_fwd_fD_Stratonovich_mean{l_ind} + data_struct.MAP_fwd_fD_Stratonovich{l_ind};
        MAP_fD_Hanggi_mean{l_ind} = MAP_fD_Hanggi_mean{l_ind} + data_struct.MAP_fD_Hanggi{l_ind};
        MAP_fD_Ito_mean{l_ind} = MAP_fD_Ito_mean{l_ind} + data_struct.MAP_fD_Ito{l_ind};
        MAP_fwd_fD_divine_mean{l_ind} = MAP_fwd_fD_divine_mean{l_ind} + data_struct.MAP_fwd_fD_divine{l_ind};
        MAP_fwd_fD_marginalized_mean{l_ind} = MAP_fwd_fD_marginalized_mean{l_ind} + data_struct.MAP_fwd_fD_marginalized{l_ind};
        fprintf('Time elapsed for this lambda: %.2f s\n', toc);
    end;
    
    
    
    
    %% Save results for this trial
    trials_data{trial} = data_struct;
 
end;

% Restore the time mesh
t_mesh = (0:N) * t_step;


%% Calculate mean
for l_ind = 1:lambda_count
    data_struct.MAP_D_mean{l_ind} = MAP_D_mean{l_ind} / trials;
%     data_struct.MAP_bck_fD_Hanggi_mean{l_ind} = MAP_bck_fD_Hanggi_mean{l_ind} / trials;
    data_struct.MAP_fD_Ito_mean{l_ind} = MAP_fD_Ito_mean{l_ind} / trials;
    data_struct.MAP_fwd_fD_Stratonovich_mean{l_ind} = MAP_fwd_fD_Stratonovich_mean{l_ind} / trials;
    data_struct.MAP_fD_Hanggi_mean{l_ind} = MAP_fD_Hanggi_mean{l_ind} / trials;
    data_struct.MAP_fwd_fD_divine_mean{l_ind} = MAP_fwd_fD_divine_mean{l_ind} / trials;
    data_struct.MAP_fwd_fD_marginalized_mean{l_ind} = MAP_fwd_fD_marginalized_mean{l_ind} / trials;
end;


% % % %% Calculate the true theoretical values in bins
% % % D_theor = cell(1, lambda_count);
% % % fD_theor = cell(1, lambda_count);
% % % for l_ind = 1:lambda_count
% % %     D_theor{l_ind} = D_func(selected_D_case, data_struct.x_bins_centers_saved{l_ind}, L);
% % %     fD_theor{l_ind} = D_func(selected_D_case, data_struct.x_bins_centers_saved{l_ind}, L)...
% % %         .* f_func(selected_f_case, data_struct.x_bins_centers_saved{l_ind}, L);
% % % end;


%% Calculate the fail rate of each prediction given the known error intervals and calculate MAP force distribution parameters
for l_ind = 1:lambda_count
    outside_count_D = zeros(1, x_bins_number_saved(l_ind));
    outside_count_divine = zeros(1, x_bins_number_saved(l_ind));
    outside_count_Ito = zeros(1, x_bins_number_saved(l_ind));
    outside_count_Stratonovich = zeros(1, x_bins_number_saved(l_ind));
    outside_count_marginalized = zeros(1, x_bins_number_saved(l_ind));
    outside_count_Hanggi = zeros(1, x_bins_number_saved(l_ind));
    combined_MAP_values_Ito = zeros(trials, x_bins_number_saved(l_ind));
    combined_MAP_values_Stratonovich = zeros(trials, x_bins_number_saved(l_ind));
    combined_MAP_values_marginalized = zeros(trials, x_bins_number_saved(l_ind));
    combined_MAP_values_divine = zeros(trials, x_bins_number_saved(l_ind));
    for trial = 1:trials
        % D
        estimate = trials_data{trial}.MAP_fwd_D{l_ind};
        outside_count_D = outside_count_D + double(estimate(4, :) > CONF_LEVEL);
%         theory = data_struct.D_theor_data{l_ind}';
%         outside_count_D = outside_count_D + double(theory < estimate(1, :) - estimate(2, :)...
%             | theory > estimate(1, :) + estimate(3, :));
        % Divine
        estimate = trials_data{trial}.MAP_fwd_fD_divine{l_ind};
        outside_count_divine = outside_count_divine + double(estimate(4, :) > CONF_LEVEL);
%         theory = data_struct.fD_theor_data{l_ind}';
%         outside_count_divine = outside_count_divine + double(theory < estimate(1, :) - estimate(2, :)...
%             | theory > estimate(1, :) + estimate(3, :));
        % Ito
        estimate = trials_data{trial}.MAP_fD_Ito{l_ind};
        outside_count_Ito = outside_count_Ito + double(estimate(4, :) > CONF_LEVEL);
%         outside_count_Ito = outside_count_Ito + double(theory < estimate(1, :) - estimate(2, :)...
%             | theory > estimate(1, :) + estimate(3, :));
        % Stratanovich
        estimate = trials_data{trial}.MAP_fwd_fD_Stratonovich{l_ind};
        outside_count_Stratonovich = outside_count_Stratonovich + double(estimate(4, :) > CONF_LEVEL);
%         outside_count_Stratonovich = outside_count_Stratonovich + double(theory < estimate(1, :) - estimate(2, :)...
%             | theory > estimate(1, :) + estimate(3, :));
        % Marginalized
        estimate = trials_data{trial}.MAP_fwd_fD_marginalized{l_ind};
        outside_count_marginalized = outside_count_marginalized + double(estimate(4, :) > CONF_LEVEL);
%         outside_count_marginalized = outside_count_marginalized + double(theory < estimate(1, :) - estimate(2, :)...
%             | theory > estimate(1, :) + estimate(3, :));
        % Hanggi
        estimate = trials_data{trial}.MAP_fD_Hanggi{l_ind};
        outside_count_Hanggi = outside_count_Hanggi + double(estimate(4, :) > CONF_LEVEL);
        %% MAP calculations
        % Divine
        combined_MAP_values_divine(trial, :) = trials_data{trial}.MAP_fwd_fD_divine{l_ind}(1, :);
        % Ito
        combined_MAP_values_Ito(trial, :) = trials_data{trial}.MAP_fD_Ito{l_ind}(1, :);
        % Stratonovich
        combined_MAP_values_Stratonovich(trial, :) = trials_data{trial}.MAP_fwd_fD_Stratonovich{l_ind}(1, :);
        % Marginalized
        combined_MAP_values_marginalized(trial, :) = trials_data{trial}.MAP_fwd_fD_marginalized{l_ind}(1, :);
    end;
    data_struct.UR_D{l_ind} = outside_count_D / trials;
    data_struct.UR_fwd_divine{l_ind} = outside_count_divine / trials;
    data_struct.UR_Ito{l_ind} = outside_count_Ito / trials;
    data_struct.UR_fwd_Stratonovich{l_ind} = outside_count_Stratonovich / trials;
    data_struct.UR_Hanggi{l_ind} = outside_count_Hanggi / trials;
    data_struct.UR_fwd_marginalized{l_ind} = outside_count_marginalized / trials;
    %% MAP calculations
    % Ito
    data_struct.MAP_fD_mean_sigma2_Ito{l_ind}(1, :) = mean(combined_MAP_values_Ito, 1);
    data_struct.MAP_fD_mean_sigma2_Ito{l_ind}(2, :) = var(combined_MAP_values_Ito, 0, 1);
    % Stratonovich
    data_struct.MAP_fD_mean_sigma2_Stratonovich{l_ind}(1, :) = mean(combined_MAP_values_Stratonovich, 1);
    data_struct.MAP_fD_mean_sigma2_Stratonovich{l_ind}(2, :) = var(combined_MAP_values_Stratonovich, 0, 1);
    % Marginalized
    data_struct.MAP_fD_mean_sigma2_marginalized{l_ind}(1, :) = mean(combined_MAP_values_marginalized, 1);
    data_struct.MAP_fD_mean_sigma2_marginalized{l_ind}(2, :) = var(combined_MAP_values_marginalized, 0, 1);
    % Divine
    data_struct.MAP_fD_mean_sigma2_divine{l_ind}(1, :) = mean(combined_MAP_values_divine, 1);
    data_struct.MAP_fD_mean_sigma2_divine{l_ind}(2, :) = var(combined_MAP_values_divine, 0, 1);

end;


%% Backup current workspace
save('backup_workspace.mat');


%% Plot
plot_article_all(data_struct);


1;























% Infer force estimates with different inference conventions


function data_struct = infer_force(data_struct, bin, trial, trials)

% Constants
load_constants;
a_ABS_MAX = 10;
bl_find_marginalized_fD_error_bars = false;		% keep off. A much faster calculation method was implemented

% Initialize
x_bins_number = data_struct.x_bins_number;
bb_prime = data_struct.MAP_bb_prime_regular_interp(bin);

[mu_n, ~, ~, ~] = get_n_parameters(bin, data_struct, 'forward');
MAP_a = zeros(x_bins_number, conventions_count, 4);


fprintf('Estimating force. Trial: %i/%i. Bin: %i/%i\n', trial, trials, bin, x_bins_number);
        % Initialize
        
%% Skip bin if empty        
bl_empty_bin = data_struct.bl_empty_bins(bin);
if bl_empty_bin
    % Save NaNs
    MAP_a(bin, enum_conv_divine, :) = ones(1, 4) * NaN;
    MAP_a(bin, enum_conv_Ito, :) = ones(1, 4) * NaN;
    MAP_a(bin, enum_conv_Stratonovich, :) = ones(1, 4) * NaN;
    MAP_a(bin, enum_conv_Hanggi, :) = ones(1, 4) * NaN;
    MAP_a(bin, enum_conv_marginalized, :) = ones(1, 4) * NaN;

    disp('Empty bin. Skipping');
end


%% Perform calculations in non-empty bins
				
%% Ito
% Prepare function
function_to_minimze = @(a) bin_a_log_posterior_func (data_struct, bin, a, 'forward');

% Make an MLE guess
MLE_guess = mu_n / t_step;

% Find confidence intervals
MAP_a(bin, enum_conv_Ito, :) = find_confidence_interval(function_to_minimze, [- a_ABS_MAX, a_ABS_MAX], true, MLE_guess,...
    CONF_LEVEL, data_struct.a_theor_data(bin), trial, bin);


%% Stratonovich (through Ito)
% Prepare function
function_to_minimze = @(a) bin_a_simple_Stratonovich_log_posterior_func(data_struct, bin, a, bb_prime, 'forward');
% Make an MLE guess
lambda = 1/2;
MLE_guess = mu_n / t_step - lambda * bb_prime;
% Find confidence intervals
MAP_a(bin, enum_conv_Stratonovich, :) = find_confidence_interval(function_to_minimze, [- a_ABS_MAX, a_ABS_MAX], true, MLE_guess,...
    CONF_LEVEL, data_struct.a_theor_data(bin), trial, bin);


%% Hanggi (through Ito)
% Prepare function
function_to_minimze = @(a) bin_a_simple_Hanggi_log_posterior_func(data_struct, bin, a, bb_prime, 'forward');
% Make an MLE guess
lambda = 1;
MLE_guess = mu_n / t_step - lambda * bb_prime;
% Find confidence intervals
MAP_a(bin, enum_conv_Hanggi, :) = find_confidence_interval(function_to_minimze, [- a_ABS_MAX, a_ABS_MAX], true, MLE_guess,...
    CONF_LEVEL, data_struct.a_theor_data(bin), trial, bin);


%% Marginalized (through Ito)
% Prepare function
function_to_minimze = @(a) bin_a_lambda_marginalized_log_posterior_func(data_struct, bin, a, bb_prime, 'forward');
% Make a guess
lambda = 0.5;
MLE_guess = mu_n / t_step - lambda * bb_prime;
% Find confidence intervals
MAP_a(bin, enum_conv_marginalized, :) = find_confidence_interval(function_to_minimze, [- a_ABS_MAX, a_ABS_MAX], bl_find_marginalized_fD_error_bars,...
    MLE_guess, CONF_LEVEL, data_struct.a_theor_data(bin), trial, bin);
        
        
% % %         %% Divine force estimate [removed in the new version of analysis code]
% % %         % Prepare function
% % %         log_function_to_minimze = @(a) bin_a_divine_inference_log_posterior_func(data_struct, trials_lambdas(trial), bin, a, bb_prime, 'forward');
% % %         % Make an MLE guess
% % %         lambda = data_struct.lambda;
% % %         MLE_guess = mu_n / t_step - lambda * bb_prime;
% % %         % Find confidence intervals if bin not empty
% % % % 		if ~
% % %         a_divine_inference = find_confidence_interval(log_function_to_minimze, [- a_ABS_MAX, a_ABS_MAX], true, MLE_guess,...
% % %             CONF_LEVEL, data_struct.a_theor_data(bin), trial, bin);
% % %         % Save
% % %         MAP_a(bin, enum_conv_divine, :) = a_divine_inference;

   
% Save results for output
data_struct.MAP_a = MAP_a;
    
    
end
    
    
    
    
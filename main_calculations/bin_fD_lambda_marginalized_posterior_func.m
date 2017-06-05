

function result = bin_fD_lambda_marginalized_posterior_func(data_struct, l_ind, bin, fD, D_grad, str_direction)

%% Globals
% global x_bins_centers_saved;


%% Constants
load_constants;
% max_D = 10;

%% Initialize
% Check grad length
if length(D_grad) > 1
    disp('Error: a scalar gradient variable required. Aborting');
    result = zeros(size(fD));
    return;
end;



%% Calculate for each value of fD independently
input_length = length(fD);
result = zeros('like', fD);
[mu_n, ~, ~, ~] = get_n_parameters(l_ind, bin, data_struct, str_direction);
% MLE_guess = (mu_n / t_step - (1/2) * D_grad) * kBT;

for ind = 1:input_length
    % Prepare the lambda limits for this specific fD
    lambda_limits = [];
    lambda_limits(1) = 0;
    count = 1;
    % The max. likely
    lambda_ML = (mu_n - fD(ind) * t_step / kBT) / (t_step * D_grad);
    if lambda_ML > 0 && lambda_ML < 1
        count = count + 1;
        lambda_limits(count) = lambda_ML;   
    end;
    count = count + 1;
    lambda_limits(count) = 1;
    % Initialize calculations
    tmp_result = 0;
    integrand_function_wrap = @(lambda) bin_mu_posterior_func(data_struct, l_ind, bin, fD(ind) * t_step / kBT + lambda * t_step * D_grad, str_direction);
    for k1 = 1:count-1
            tmp_result = tmp_result + integral(integrand_function_wrap, lambda_limits(k1), lambda_limits(k1 + 1),...
                'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE);
    end;
    result(ind) = tmp_result;
end;


% Rescale
result = result * t_step / kBT;











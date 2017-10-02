

function result = bin_a_lambda_marginalized_posterior_func(data_struct, bin, a, bb_prime, str_direction)

%% Globals
% global x_bins_centers_saved;


%% Constants
load_constants;
% max_D = 10;

%% Initialize
% Check grad length
if length(bb_prime) > 1
    disp('Error: a scalar gradient variable required. Aborting');
    result = zeros(size(a));
    return;
end;



%% Calculate for each value of fD independently
input_length = length(a);
result = zeros('like', a );
[mu_n, ~, ~, ~] = get_n_parameters(bin, data_struct, str_direction);
% MLE_guess = (mu_n / t_step - (1/2) * D_grad) * kBT;

for ind = 1:input_length
    % Prepare the lambda limits for this specific a: if the most likely lambda is found within the [0;1] interval, 
	% split the integral in 2 takin the ML lambda into account
    lambda_limits = [];
    lambda_limits(1) = 0;
    count = 1;
    % Most likely lambda
    lambda_ML = (mu_n /t_step - a(ind)) / bb_prime;
	% Next check automatically takes into account 0 gradient
	if lambda_ML > 0 && lambda_ML < 1
		count = count + 1;
        lambda_limits(count) = lambda_ML;   
    end;
	count = count + 1;
    lambda_limits(count) = 1;
    % Initialize calculations
    tmp_result = 0;
    integrand_function_wrap = @(lambda) bin_mu_posterior_func(data_struct, bin, t_step .* (a(ind) + lambda * bb_prime), str_direction);
    for k1 = 1:count-1
            tmp_result = tmp_result + integral(integrand_function_wrap, lambda_limits(k1), lambda_limits(k1 + 1),...
                'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE);
    end;
    result(ind) = tmp_result;
end;


% Rescale
result = result * t_step;











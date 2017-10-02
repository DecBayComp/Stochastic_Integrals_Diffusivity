

function result = bin_a_divine_inference_log_posterior_func(data_struct, lambda, bin, a, bb_prime, str_direction)

%% Globals


%% Constants
load_constants;


%% Initialize
% Check grad length
if length(bb_prime) > 1
    error('A scalar gradient variable required. Aborting');
    result = zeros(size(a));
    return;
end;


%% Calculate
result = log(t_step) + bin_mu_log_posterior_func(data_struct, bin, t_step .* (a + lambda * bb_prime), str_direction);










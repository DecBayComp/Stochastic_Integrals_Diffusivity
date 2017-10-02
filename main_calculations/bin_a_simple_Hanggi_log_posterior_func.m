

function result = bin_a_simple_Hanggi_log_posterior_func(data_struct, bin, a, bb_prime, str_direction)


%% Constants
load_constants;


%% Initialize
% Check grad length
if length(bb_prime) > 1
    disp('Error: a scalar gradient variable required. Aborting');
    result = zeros(size(a));
    return;
end;


%% Calculate
lambda = 1;
result = log(t_step) + bin_mu_log_posterior_func(data_struct, bin, t_step .* (a + lambda * bb_prime), str_direction);









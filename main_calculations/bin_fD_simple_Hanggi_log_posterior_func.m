

function result = bin_fD_simple_Stratonovich_log_posterior_func(data_struct, bin, fD, D_grad, str_direction)


%% Constants
load_constants;


%% Initialize
% Check grad length
if length(D_grad) > 1
    disp('Error: a scalar gradient variable required. Aborting');
    result = zeros(size(fD));
    return;
end;


%% Calculate
lambda = 1;
result = log(t_step / kBT) + bin_mu_log_posterior_func(data_struct, bin, fD .* t_step ./ kBT + lambda * t_step * D_grad, str_direction);









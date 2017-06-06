

function result = bin_fD_divine_inference_log_posterior_func(data_struct, lambda, bin, fD, D_grad, str_direction)

%% Globals


%% Constants
load_constants;


%% Initialize
% Check grad length
if length(D_grad) > 1
    error('A scalar gradient variable required. Aborting');
end;


%% Calculate
result = log(t_step / kBT) + bin_mu_log_posterior_func(data_struct, bin, fD .* t_step ./ kBT + lambda * t_step * D_grad, str_direction);










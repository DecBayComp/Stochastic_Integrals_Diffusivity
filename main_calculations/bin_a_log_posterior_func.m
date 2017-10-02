

%% Defining the sigma distribution function for each individual bin
function out =  bin_a_log_posterior_func (data_struct, bin, fD, str_direction)


%% Constants
load_constants;


%% Calculating
out = log(t_step / kBT) + bin_mu_log_posterior_func (data_struct, bin, fD .* t_step ./ kBT, str_direction);








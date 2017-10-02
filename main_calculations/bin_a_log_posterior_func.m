

%% Defining the sigma distribution function for each individual bin
function out =  bin_a_log_posterior_func (data_struct, bin, a, str_direction)


%% Constants
load_constants;


%% Calculating
out = log(t_step) + bin_mu_log_posterior_func (data_struct, bin, t_step .* (a + 0), str_direction);








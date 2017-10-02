

%% Defining the sigma distribution function for each individual bin
function out =  bin_a_posterior_func (data_struct, bin, a, str_direction)

%% Globals

%% Calculating
out = exp(bin_a_log_posterior_func(data_struct, bin, a, str_direction));







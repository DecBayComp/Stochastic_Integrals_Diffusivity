

%% Defining the sigma distribution function for each individual bin
function out =  bin_fD_posterior_func (data_struct, bin, fD, str_direction)

%% Globals

%% Calculating
out = exp(bin_fD_log_posterior_func(data_struct, bin, fD, str_direction));







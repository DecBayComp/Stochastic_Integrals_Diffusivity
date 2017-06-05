

%% Defining the sigma distribution function for each individual bin
function out =  bin_mu_posterior_func (data_struct, l_ind, bin, mu, str_direction)

out = bin_mu_log_posterior_func (data_struct, l_ind, bin, mu, str_direction);

out = exp(out);




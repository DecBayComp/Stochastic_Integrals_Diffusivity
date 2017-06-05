

function out =  bin_sigma_squared_posterior_func (l_ind, bin, sigma2, data_struct_n, str_direction)


out = exp(bin_sigma_squared_log_posterior_func (l_ind, bin, sigma2, data_struct_n, str_direction));




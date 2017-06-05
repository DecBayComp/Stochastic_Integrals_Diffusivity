

function out =  bin_D_posterior_func (l_ind, bin, D, t_step, data_struct_n, str_direction)


%% Constants
% load_constants;


%% Calculate parameters
out = exp(bin_D_log_posterior_func (l_ind, bin, D, t_step, data_struct_n, str_direction));


1;
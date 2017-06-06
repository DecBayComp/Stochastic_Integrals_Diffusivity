

function out =  bin_D_log_posterior_func (bin, D, t_step, data_struct, str_direction)


%% Constants
% load_constants;


%% Calculate parameters
out = log(2 * t_step) + bin_sigma_squared_log_posterior_func (bin, 2 * t_step * D, data_struct, str_direction);


1;
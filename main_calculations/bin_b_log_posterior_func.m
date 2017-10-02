

function out =  bin_b_log_posterior_func (bin, b, t_step, data_struct, str_direction)


%% Constants
% load_constants;


%% Calculate parameters
out = log(2 * b * t_step) + bin_sigma_squared_log_posterior_func (bin, b^2 * t_step, data_struct, str_direction);


1;
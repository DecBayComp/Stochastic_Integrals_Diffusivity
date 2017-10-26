%% Calculate log probability density of b=sqrt(2*D)
% Checked: OK


function out =  bin_b_posterior_func (bin, b, t_step, data_struct, str_direction)


%% Calculate parameters
out = exp(bin_b_log_posterior_func (bin, b, t_step, data_struct, str_direction));

% Set the posterior to 0 for b<0 (even if it is already set by the log_b posterior function)
out(b<=0) = 0;


1;


function out =  bin_sigma_squared_log_posterior_func (bin, sigma2, data_struct_n, str_direction)


%% Constants
% load_constants;

% % Check for negative input values
% if 


%% Calculate parameters
[~, ~, nu_n, sigma2_n] = get_n_parameters(bin, data_struct_n, str_direction);


% Calculate log of the output
out = nu_n/2 * log(nu_n .* sigma2_n ./ 2) - gammaln(nu_n/2) - (nu_n/2 + 1) .* log(sigma2)...
    - nu_n .* sigma2_n ./ (2 * sigma2);

% Generate -Inf output for input non-posititve values of sigma2
% D = 0 is also unphysical
out(sigma2 <= 0) = -Inf;


1;









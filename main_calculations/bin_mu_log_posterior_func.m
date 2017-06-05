

%% Defining the sigma distribution function for each individual bin
function out =  bin_mu_log_posterior_func (data_struct, l_ind, bin, mu, str_direction)

%% Globals
% global kBT;
% global t_step;

%% Calculate parameters
% _n parameters
[mu_n, kappa_n, nu_n, sigma2_n] = get_n_parameters(l_ind, bin, data_struct, str_direction);

%% Calculating
out = 1/2 * log(kappa_n / (pi * nu_n * sigma2_n))...
    + gammaln(nu_n/2) - gammaln((nu_n - 1)/2)...
    - nu_n/2 * log(1 + kappa_n * (mu - mu_n).^2 / (nu_n * sigma2_n));





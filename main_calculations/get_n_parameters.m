%% Calculate '_n' parameters of the normal-inverse chi-squared function to perform inference


function [mu_n, kappa_n, nu_n, sigma2_n] = get_n_parameters(bin, data_struct_n, str_direction)



%% Load data
V_j = data_struct_n.V_j;
n = data_struct_n.n_j(bin);
dx_mean_bin = data_struct_n.dx_mean_in_bins(bin);



%% Calculate '_pi' parameters
[mu_pi, kappa_pi, nu_pi, sigma2_pi] = get_pi_parameters(data_struct_n);



%% Calculate '_n' parameters
mu_n = (n * dx_mean_bin + kappa_pi * mu_pi) / (n + kappa_pi);
kappa_n = n + kappa_pi;
nu_n = n + nu_pi;
sigma2_n = (n * kappa_pi * (mu_pi - dx_mean_bin)^2/(n + kappa_pi) + n * V_j(bin) ...
    + nu_pi * sigma2_pi) / (n + nu_pi);

1;




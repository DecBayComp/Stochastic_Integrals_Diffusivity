%% Calculate '_c' parameters of the normal-inverse chi-squared function for calculation of the Bayes factor
% I will need to make sure that _pi parameters are the same ones as what I use for the other prior!


function [mu_c, kappa_c, nu_c, sigma2_c] = get_c_parameters(bin, data_struct)



%% Load data
V_j = data_struct.V_j(bin);
n = data_struct.n_j(bin);
dx_mean_bin = data_struct.dx_mean_in_bins(bin);



%% Calculate '_pi' parameters
[~, ~, nu_pi, sigma2_pi] = get_pi_parameters(data_struct);



%% Calculate '_n' parameters
mu_c = dx_mean_bin;
kappa_c = n;
nu_c = n + nu_pi;
sigma2_c = (n * V_j + nu_pi * sigma2_pi) / (n + nu_pi);

1;




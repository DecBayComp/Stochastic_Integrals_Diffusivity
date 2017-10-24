%% Calculate '_n' parameters of the normal-inverse chi-squared function used as prior



function [mu_pi, kappa_pi, nu_pi, sigma2_pi] = get_pi_parameters(data_struct)



%% Load data
dx_Mean = data_struct.dx_mean_all_bins_all_trials;
V = data_struct.V;


% _pi parameters
mu_pi = dx_Mean;
kappa_pi = 1;
nu_pi = 1e0;
sigma2_pi = (2 + nu_pi) / nu_pi * V;

1;




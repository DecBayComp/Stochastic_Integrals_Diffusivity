function [mu_n, kappa_n, nu_n, sigma2_n] = get_n_parameters(l_ind, bin, data_struct_n, str_direction)

%% Constants and globals
% load_constants;
% global dx_mean_in_bins_saved;
% global dx_Mean;
% global n_j;
% global V_j;
% global V;


%% Initialize
% Identify the direction
if strcmp(str_direction, 'forward')
    V_j = data_struct_n.V_fwd_j;
    dx_mean_in_bins_saved = data_struct_n.dx_fwd_mean_in_bins_saved;
elseif strcmp(str_direction, 'backward')
    V_j = data_struct_n.V_bck_j;
    dx_mean_in_bins_saved = data_struct_n.dx_bck_mean_in_bins_saved;
else
    error('The direction should either be ''forward'' or ''backward''');
end;
%
n_j = data_struct_n.n_j;
n = n_j{l_ind}(bin);
dx_Mean_bin = dx_mean_in_bins_saved{l_ind}(bin);

% _pi parameters
% % % data_struct_pi = struct('dx_Mean', dx_Mean, 'V', V);
[mu_pi, kappa_pi, nu_pi, sigma2_pi] = get_pi_parameters(l_ind, data_struct_n);

% _n parameters
mu_n = (n * dx_Mean_bin + kappa_pi * mu_pi) / (n + kappa_pi);
kappa_n = n + kappa_pi;
nu_n = n + nu_pi;
sigma2_n = (n * kappa_pi * (mu_pi - dx_Mean_bin)^2/(n + kappa_pi) + n * V_j{l_ind}(bin) ...
    + nu_pi * sigma2_pi) / (n + nu_pi);

1;




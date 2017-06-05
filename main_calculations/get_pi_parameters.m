function [mu_pi, kappa_pi, nu_pi, sigma2_pi] = get_pi_parameters(l_ind, data_struct)


% %% Constants and globals
% if length(varargin) > 0
% else
%     global dx_Mean;
%     global V;
% end;

dx_Mean = data_struct.dx_Mean;
V = data_struct.V;


% _pi parameters
mu_pi = dx_Mean(l_ind);
kappa_pi = 1;
nu_pi = 1e0;
sigma2_pi = (2 + nu_pi) / nu_pi * V(l_ind);

1;




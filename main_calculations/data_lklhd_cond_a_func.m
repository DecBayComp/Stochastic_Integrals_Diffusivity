%% Calculate data likelihood for the model iwth different diffusivity and diffusivity gradients across bins, but same force.
% The log of the a-independent pre-factor is given independently to simplify integration
% The initial parameters (n, V, grad, mu_c, kappa_c, nu_c, sigma2_c) should be vectors of the same length
% Treat a as a vector of values.
% a-related output is vectors. Pre-factors are scalars
%
% Output: 
% a_part - the product of a-dependent part across all bins
% a_prior - a-dependent part of the a prior
% ln_pre_factor - a-independent part of the product across bins
% ln_pre_factor_a_prior - a-independent part of the a prior (dt / sqrt(2*pi*V))




function [a_part, a_prior, log_pre_factor, log_pre_factor_a_prior] = data_lklhd_cond_a_func(a, n, V, grad, mu_c, kappa_c, nu_c, sigma2_c, nu_pi, sigma2_pi, lambda, t_step, dx_mean_tot, V_tot)


% % % %% Constants
% % % nu_pi = 1;
% % % sigma

% Calculate bin number
m = length(n);

%% Calculate the contribution of data likelihood
a_part = ones(size(a));
log_pre_factor = 1;
for j=1:m
	% Use different expressions if the bin has zero points, one point ore more points
	n_j = n(j);
	if n_j == 0
		continue;
	elseif n_j == 1
		a_part = a_part .* (1 + ((a + lambda * grad(j)) * t_step - mu_c(j)).^2 / nu_pi / sigma2_pi) .^ (-(nu_pi + 1) / 2);
		
		log_pre_factor = log_pre_factor + gammaln((1+ nu_pi)/2) - 1/2 * log(pi * nu_pi * sigma2_pi) - gammaln(nu_pi / 2);
	elseif n_j > 1
		a_part = a_part .* (1 + kappa_c(j) * ((a + lambda * grad(j)) * t_step - mu_c(j)).^2 / nu_c(j) / sigma2_c(j)) .^ (-(nu_c(j) + 1) / 2);
		
		log_pre_factor = log_pre_factor + gammaln((n_j - 3)/2) + gammaln((1 + nu_c(j)) / 2) - gammaln(nu_c(j) / 2) + 1/2 * log(kappa_c(j))...
			- log(2) - 1/2 * log(n_j) - (n_j-2)/2 * log(pi) - (n_j-3)/2 * log(n_j * V(j)) - 1/2 * log(pi * nu_c(j) * sigma2_c(j));
	end;
end;



%% Calculate the part from the a prior
a_prior = exp(-(a * t_step - dx_mean_tot) .^ 2 / 2 / V_tot);
log_pre_factor_a_prior = log(t_step) - 1/2 * log(2 * pi * V_tot);

1;






















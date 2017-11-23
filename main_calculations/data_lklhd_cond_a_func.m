%% Calculate data likelihood for the model with different diffusivity and diffusivity gradients across bins, but same force
% The log of the a-independent pre-factor in bins and in prior is given independently to simplify integration
% The initial parameters (n, V, grad, mu_c, kappa_c, nu_c, sigma2_c) should be vectors of the same length
% a is treated as a vector of values of arbitrary length
% a-related output is vectors. Pre-factors are scalars
%
% Output: 
% a_lklhd_bins - the product of a-dependent parts of likelihood across all bins
% a_prior - a-dependent part of the a prior
% log_pre_factor_bins - a-independent part of the likelihood product across bins
% log_pre_factor_a_prior - a-independent part of the a prior



function [a_lklhd_bins, a_prior, log_pre_factor_bins, log_pre_factor_prior] = data_lklhd_cond_a_func(a, n, V, grad, mu_c, kappa_c, nu_c, sigma2_c,...
																		nu_pi, sigma2_pi, lambda, t_step, dx_mean_tot, V_tot)

%% Initialize
m = length(n);
a_lklhd_bins = ones(size(a));



%% Calculate unnormalized a likelihood in bins and the pre-factor from bin likelihood
log_pre_factor_bins = 1;
for j=1:m
	% Use different expressions if the bin has zero points, one point ore more points
	n_j = n(j);
	if n_j == 0
		continue;
		
	elseif n_j == 1
		a_lklhd_bins = a_lklhd_bins .* (1 + ((a + lambda * grad(j)) * t_step - mu_c(j)).^2 / nu_pi / sigma2_pi) .^ (-(nu_pi + 1) / 2);
		
		log_pre_factor_bins = log_pre_factor_bins + gammaln((1+ nu_pi)/2) - 1/2 * log(pi * nu_pi * sigma2_pi) - gammaln(nu_pi / 2);
		
	elseif n_j > 1
		a_lklhd_bins = a_lklhd_bins .* (1 + kappa_c(j) * ((a + lambda * grad(j)) * t_step - mu_c(j)).^2 / nu_c(j) / sigma2_c(j)) .^ (-(nu_c(j) + 1) / 2);
		
		log_pre_factor_bins = log_pre_factor_bins + gammaln((n_j - 3)/2) + gammaln((1 + nu_c(j)) / 2) - gammaln(nu_c(j) / 2) + 1/2 * log(kappa_c(j))...
			- log(2) - 1/2 * log(n_j) - (n_j-2)/2 * log(pi) - (n_j-3)/2 * log(n_j * V(j)) - 1/2 * log(pi * nu_c(j) * sigma2_c(j));
	end;
end;



%% Calculate the a prior parts
a_prior = exp(-(a * t_step - dx_mean_tot) .^ 2 / 2 / V_tot);
log_pre_factor_prior = log(t_step) - 1/2 * log(2 * pi * V_tot);






















%% Calculate Bayes factors for the marginalized and fixed-lambda models with force relative to corresponding no-force models
% The function uses the MAP value of the diffusivity gradient bb' in the bin.
% K is > 1 if there is a local force, i.e. K = Pr(force_model) / Pr(no_force_model)
%
% Output: ??


function log_K = calculate_bayes_factor_new(data_struct)

%% Constants
load_constants;
ABS_TOLERANCE = 1e-7;
REL_TOLERANCE = 1e-7;



%% Initialize
lambda_true = data_struct.lambda;
log_K = zeros(conventions_count, 1);
lambda_array = [0, 0.5, 1, lambda_true];

% Load mean jump and variance across all bins
dx_mean_tot = data_struct.dx_mean_all_bins;
V_tot = data_struct.V;

% Load individual bin parameters and reshape into vectors
n = reshape(data_struct.n_j, [], 1);
V = reshape(data_struct.V_j, [], 1);
grad = reshape(data_struct.MAP_bb_prime_regular_interp, [], 1);

% _pi parameters
nu_pi = 1;
sigma2_pi = V_tot;

% _c parameters for the posteriors in each bin
mu_c = reshape(data_struct.dx_mean_in_bins, [], 1);
kappa_c = n;
nu_c = n + nu_pi;
sigma2_c = (n .* V + nu_pi * sigma2_pi) ./ (n + nu_pi);



%% Calculate the Bayes factor for fixed-lambda conventions
for convention = 1:conventions_count-1
	% Set inference lambda
	lambda = lambda_array(convention);

	% Wrap likelihood function to return only a including both bin likelihood and a prior
	data_lklhd_cond_a_func_wrap = @(a) data_lklhd_cond_a_func_a_integrand(a, n, V, grad, mu_c, kappa_c, nu_c, sigma2_c, nu_pi, sigma2_pi,...
		lambda, t_step, dx_mean_tot, V_tot);

	% Separate integration region [-Inf, Inf] into separate zones by adding the most likely a values in each bin
	a_ML = - lambda * grad;
	a_limits = sort([a_ML; -Inf; Inf]);

	% To calculate the probability of the force model, marginalize the unnormalized posterior over a 
	a_integral = 0;
	for i = 1:length(a_limits) - 1
		a_integral = a_integral + integral(data_lklhd_cond_a_func_wrap, a_limits(i), a_limits(i+1), 'RelTol', REL_TOLERANCE, 'AbsTol', ABS_TOLERANCE);
	end;

	% For probability of the null-model, take the likelihood part without the 'a' prior
	[a0_part, ~, ~, log_pre_factor_a_prior] = data_lklhd_cond_a_func(0, n, V, grad, mu_c, kappa_c, nu_c, sigma2_c, nu_pi, sigma2_pi,...
		lambda, t_step, dx_mean_tot, V_tot);

	% Combine into unnormalized probabilities of the force and no-force models
	log_force_model_unnorm_probability = log(a_integral) + log_pre_factor_a_prior;
	log_null_model_unnorm_probability = log(a0_part);
	
	% Calculate the Bayes factor as the ratio of the model probabilites. The force model is upstairs
	log_K(convention) = log_force_model_unnorm_probability - log_null_model_unnorm_probability;
end;









% % % for bin = 1:x_bins_number
% % % 	% Calculate the combined ('_c') N-inv-chi2 parameters for the bin
% % % 	[mu_c, kappa_c, nu_c, sigma2_c] = get_c_parameters(bin, data_struct);
% % % 
% % % 	% Load the MAP value of the gradient
% % % 	MAP_bb_prime = data_struct.MAP_bb_prime_regular_interp(bin);
% % % 	% Determine gradient sign
% % % 	sgn_bb_prime = sign(MAP_bb_prime);
% % % 
% % % 	
% % % 	%% Fixed-lambda estimates
% % % 	% Pre-factor
% % % 	ln_pre_factor = 1/2 * log(pi * nu_c * sigma2_c / kappa_c) + gammaln(nu_c/2) - gammaln((1+nu_c)/2);
% % % 	
% % % 	
% % % 	%% Non-marginalized models
% % % 	% Lambda-factor
% % % 	lambda = [0, 0.5, 1, lambda_true];
% % % 	z = (lambda * MAP_bb_prime * t_step - mu_c) * sqrt(kappa_c / nu_c / sigma2_c);
% % % 	ln_lambda_factor = (nu_c + 1)/2 * log(1 + z.^2);
% % % 	
% % % 	% Save result
% % % 	log_K([enum_conv_Ito, enum_conv_Stratonovich, enum_conv_Hanggi, enum_conv_divine], bin) = ln_pre_factor + ln_lambda_factor;
% % % 		
% % % 	
% % % 	
% % % 	%% Marginalized model
% % % 	% Pre-factor of the Bayes factor
% % % 	ln_pre_factor = gammaln(nu_c/2) - gammaln((1+nu_c)/2) + log (MAP_bb_prime * t_step * sqrt(pi) * sgn_bb_prime);
% % % 	
% % % 	% Integration limits for the 2F1 hypergeometric function
% % % 	z1 = - mu_c * sqrt(kappa_c / nu_c / sigma2_c);
% % % 	z2 = (MAP_bb_prime * t_step - mu_c) * sqrt(kappa_c / nu_c / sigma2_c);
% % % 	
% % % 	% Add 0 to integration points if it lies between the limits
% % % 	if z1 * z2 < 0
% % % 		z_int_points = [z1, 0, z2];
% % % 	else
% % % 		z_int_points = [z1, z2];
% % % 	end;
% % % 	
% % % 	% Prepare the function to integrate over lambda
% % % 	integrand_func = @(z) (1+z.^2).^(-(nu_c + 1)/2);
% % % 
% % % 	% Split integral at 0 to increase accuracy
% % % 	hyperg = 0;
% % % 	for i = 1:length(z_int_points) - 1
% % % 		hyperg = hyperg + integral(integrand_func, z_int_points(i), z_int_points(i+1), 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE);
% % % 	end;
% % % 	ln_hyperg = log(hyperg * sgn_bb_prime);
% % % 
% % % 	% Combine both into the Bayes factor for the marginalized model
% % % 	log_K(enum_conv_marginalized, bin) = ln_pre_factor - ln_hyperg;
% % % end;



















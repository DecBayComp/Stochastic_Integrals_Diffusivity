%% Calculate Bayes factors for the marginalized and fixed-lambda models with force relative to corresponding no-force models
% The function uses the MAP value of the diffusivity gradient bb' in the bin.
% K is > 1 if there is a local force, i.e. K = Pr(force_model) / Pr(no_force_model)
%
% Output: ??


function [log_K_L, log_K_G] = calculate_bayes_factor(data_struct)

%% Constants
load_constants;
ABS_TOLERANCE = 1e-7;
REL_TOLERANCE = 1e-7;



%% Initialize
lambda_true = data_struct.lambda;
lambda_array = [0, 0.5, 1, lambda_true];

% Load mean jump and variance across all bins
dx_mean_tot = data_struct.dx_mean_all_bins;
V_all = data_struct.V_all_bins;

% Load individual bin parameters and reshape into vectors
n = reshape(data_struct.n_j, [], 1);
dx_mean = reshape(data_struct.dx_mean_in_bins, [], 1); 
V = reshape(data_struct.V_j, [], 1);
grad = reshape(data_struct.MAP_bb_prime_regular_interp, [], 1);

x_bins_count = length(n);
m = x_bins_count;

% Determine gradient sign
sgn_grad = sign(grad);

% _pi parameters
n_pi = 4;
mu_pi = dx_mean_tot;
sigma2_pi = V_all;

% _0 and _a parameters independent of lambda 
n_0 = n + n_pi;
n_a = n_0;

mu_a = (n .* dx_mean + n_pi .* mu_pi) ./ n_a;
sigma2_a = (n .* dx_mean.^2 + n_pi * mu_pi.^2 + n .* V + n_pi .* sigma2_pi) ./ n_a - mu_a.^2;
sigma2_0_func = @(lambda) (n .* (lambda * grad * t_step - dx_mean) .^ 2 + n .* V + n_pi * sigma2_pi) ./ n_0;

log_A_0 = log(2) + n_pi/2 * log(2*pi) - log(n_pi * sigma2_pi) - gammaln(n_pi/2 - 1) + n_pi /2 * log(n_pi * sigma2_pi / 2);
log_A_a = log (2 * pi) + (n_pi - 3)/2 * log(pi * n_pi * sigma2_pi) - gammaln((n_pi - 3)/2);

% Integration limits for the 2F1 hypergeometric function
y1 = - sqrt(n ./ (n .* V + n_pi * sigma2_pi)) .* dx_mean;
y2 =  sqrt(n ./ (n .* V + n_pi * sigma2_pi)) .* (grad * t_step - dx_mean);

% Bayes pre-factors for the local and global marginalized model
log_K_L_marg_prefactor = log_A_a - log_A_0 + gammaln((n_a - 3)/2) - gammaln(n_0/2 - 1) + (n_0 - n_a + 1)/2 * log(pi) - 1/2 * log(n_a) ...
		+ (n_0 / 2 - 1) .* log(n .* V + n_pi * sigma2_pi) - (n_a - 3)/2 .* log(n_a .* sigma2_a);

tmp =  gammaln((n_a - 3) / 2) - gammaln(n_0 / 2 - 1) + (n_0 / 2 - 1) .* log(n .* V + n_pi * sigma2_pi) - (n_a - 3) / 2 .* log (n_a .* sigma2_a) ...
	+ (n_0 - n_a + 1) / 2 * log(pi) - 1/2 * log(n_a);
log_K_G_marg_prefactor = m * (log_A_a - log_A_0) + sum(tmp);

	
% Hypergeometric function pre-factor for the marginalized model
log_marg_integral_prefactor = 1/2 * log((n .* V + n_pi * sigma2_pi) ./ n) - log(t_step * grad .* sgn_grad);



%% Calculate the local and global Bayes factors for fixed-lambda conventions
log_K_L = zeros(conventions_count, x_bins_count);
log_K_G = zeros(conventions_count, 1);
for convention = 1:conventions_count-1
	% Set inference lambda
	lambda = lambda_array(convention);
	
	% Calculate sigma2_0 for the given lambda
	sigma2_0 = sigma2_0_func(lambda);
	
	% Calculate local Bayes factor for fixed-lambda models
	log_K_L(convention, :) = log_A_a - log_A_0 + gammaln((n_a - 3)/2) - gammaln(n_0/2 - 1) + (n_0 - n_a + 1)/2 * log(pi) - 1/2 * log(n_a) ...
		+ (n_0 / 2 - 1) .* log(n_0 .* sigma2_0) - (n_a - 3)/2 .* log(n_a .* sigma2_a);
	
	
	
	%% Calculate local Bayes factor for the marginalized model in each bin
	for bin = 1:x_bins_count
		% Add 0 to integration points list if it lies between the limits
		if y1(bin) * y2(bin) < 0
			int_limits = [y1(bin), 0, y2(bin)];
		else
			int_limits = [y1(bin), y2(bin)];
		end;

		% Prepare the function to integrate over lambda in the current bin
		lambda_integrand_func = @(y) (1 + y.^2) .^ (1 - n_0(bin)/2);

		% Split integral at 0 to increase accuracy and calculate
		hyperg = 0;
		for i = 1:length(int_limits) - 1
			hyperg = hyperg + integral(lambda_integrand_func, int_limits(i), int_limits(i+1), 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE);
		end;
		log_hyperg = log(hyperg * sgn_grad(bin));
		
		log_K_L(enum_conv_marginalized, bin) = log_K_L_marg_prefactor(bin) + log_marg_integral_prefactor (bin) - log_hyperg;
	end;
	
	
	
% % % 	%% Calculate the global Bayes factor for fixed-lambda models
% % % 	% Prepare the array of values in individual bins
% % % 	tmp =  gammaln((n_a - 3) / 2) - gammaln(n_0 / 2 - 1) + (n_0 / 2 - 1) .* log(n_0 .* sigma2_0) - (n_a - 3) / 2 .* log (n_a .* sigma2_a) + (n_0 - n_a + 1) / 2 * log(pi) - 1/2 * log(n_a);
% % % 	log_K_G(convention) = m * log_A_a - m * log_A_0 + sum(tmp);
% % % 	
	
	
	
	
	

% % % % % 	% Wrap likelihood function to return only a including both bin likelihood and a prior
% % % % % 	data_lklhd_cond_a_func_wrap = @(a) data_lklhd_cond_a_func_a_integrand(a, n, V, grad, mu_c, kappa_c, nu_c, sigma2_c, nu_pi, sigma2_pi,...
% % % % % 		lambda, t_step, dx_mean_tot, V_all);
% % % % % 
% % % % % 	% Separate integration region [-Inf, Inf] into separate zones by adding the most likely a values in each bin
% % % % % 	a_ML = - lambda * grad;
% % % % % 	a_limits = sort([a_ML; -Inf; Inf]);
% % % % % 
% % % % % 	% To calculate the probability of the force model, marginalize the unnormalized posterior over a 
% % % % % 	a_integral = 0;
% % % % % 	for i = 1:length(a_limits) - 1
% % % % % 		a_integral = a_integral + integral(data_lklhd_cond_a_func_wrap, a_limits(i), a_limits(i+1), 'RelTol', REL_TOLERANCE, 'AbsTol', ABS_TOLERANCE);
% % % % % 	end;
% % % % % 
% % % % % 	% For probability of the null-model, take the likelihood part without the 'a' prior
% % % % % 	[a0_part, ~, ~, log_pre_factor_a_prior] = data_lklhd_cond_a_func(0, n, V, grad, mu_c, kappa_c, nu_c, sigma2_c, nu_pi, sigma2_pi,...
% % % % % 		lambda, t_step, dx_mean_tot, V_all);
% % % % % 
% % % % % 	% Combine into unnormalized probabilities of the force and no-force models
% % % % % 	log_force_model_unnorm_probability = log(a_integral) + log_pre_factor_a_prior;
% % % % % 	log_null_model_unnorm_probability = log(a0_part);
% % % % % 	
% % % % % 	% Calculate the Bayes factor as the ratio of the model probabilites. The force model is upstairs
% % % % % 	log_K_L(convention) = log_force_model_unnorm_probability - log_null_model_unnorm_probability;
end;

%% Calculate the global Bayes factor for fixed-lambda models by summing the local Bayes factor over bins
log_K_G(1:end-1) = sum(log_K_L(1:end-1, :), 2);



%% Calculate the global Bayes factor integral for the marginalized model
% Get handle to the function to integrate over lambda
lambda_integrand_func = get_global_marg_integrand_hand(n, grad, dx_mean, V, n_pi, sigma2_pi, n_0, t_step);

% Separate integration region lambda \in [0, 1] into separate zones by adding the most likely a values in each bin
lambda_ML = dx_mean ./ grad / t_step;

% Keep only those that are within the region
lambda_ML = lambda_ML(lambda_ML >= 0 & lambda_ML <= 1);

% Combine with original borders
lambda_limits = sort([lambda_ML; 0; 1]);

% Integrate
lambda_integral = 0;
for i = 1:length(lambda_limits) - 1
	lambda_integral = lambda_integral + integral(lambda_integrand_func, lambda_limits(i), lambda_limits(i+1), 'RelTol', REL_TOLERANCE, 'AbsTol', ABS_TOLERANCE);
end;



%% Combine the lambda integral with the pre-factor for the global Bayes factor of the marginalized model
log_K_G(enum_conv_marginalized) = log_K_G_marg_prefactor - log(lambda_integral);
	

% % % disp (log_K_G);

% disp (sum(log_K_L, 2));

end



%% Define the function to integrate over lambda for the global marginalized model
function hand = get_global_marg_integrand_hand(n, grad, dx_mean, V, n_pi, sigma2_pi, n_0, t_step)
	% Create handle to the integrand function as a function of only lambda
	hand = @ integrand_func;

	% Define the integrand function. It shares variables with the envelopping function
	function vals = integrand_func(lambda_array)
		% Initialize
		vals = zeros(size(lambda_array));
		lambda_count = length(lambda_array);

		% Parse lambda values supplied to the function
		for lambda_ind = 1:lambda_count
			lambda = lambda_array(lambda_ind);
			
			% Calculate logarithms of individual bins contributions
			log_vals = (1 - n_0 / 2) .* log(1 + n .* (lambda .* grad * t_step - dx_mean) .^ 2 ./ (n .* V + n_pi .* sigma2_pi));
			
			% Calculate product
			vals(lambda_ind) = exp(sum(log_vals));
		end
	end
end







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



















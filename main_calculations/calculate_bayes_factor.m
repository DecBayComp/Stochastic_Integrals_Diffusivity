%% Calculate Bayes factor for a given bin to compare the marginalized model with no force and the marginalized model with a local force.
% Both include the spurious force.
% The function uses the MAP value of the bb' in the bin.
% The no-force model in downstairs in the ratio
% So in total K is larger if there is local force, i.e. K = Pr(force_model) / Pr(no_force_model)


function K = calculate_bayes_factor(data_struct)



%% Constants
load_constants;



%% Initialize
x_bins_centers = data_struct.x_bins_centers;
x_bins_number = length(data_struct.x_bins_centers);
log_K = zeros(conventions_count, x_bins_number);
lambda_true = data_struct.lambda;


for bin = 1:x_bins_number
	% Calculate the combined ('_c') N-inv-chi2 parameters
	[mu_c, kappa_c, nu_c, sigma2_c] = get_c_parameters(bin, data_struct);

	% Load the MAP value of the gradient
	MAP_bb_prime = data_struct.MAP_bb_prime_regular_interp(bin);
	% Determine gradient sign
	sgn_bb_prime = sign(MAP_bb_prime);

	
	%% Fixed-lambda estimates
	% Pre-factor
	ln_pre_factor = 1/2 * log(pi * nu_c * sigma2_c / kappa_c) + gammaln(nu_c/2) - gammaln((1+nu_c)/2);
	
	
	%% Ito
	% Lambda-factor
	lambda = [0, 0.5, 1, lambda_true];
	z = (lambda * MAP_bb_prime * t_step - mu_c) * sqrt(kappa_c / nu_c / sigma2_c);
	ln_lambda_factor = (nu_c + 1)/2 * log(1 + z.^2);
	log_K([enum_conv_Ito, enum_conv_Stratonovich, enum_conv_Hanggi, enum_conv_divine], bin) = ln_pre_factor + ln_lambda_factor;
	
	
	
	
	%% Marginalized integral calculation
	% Pre-factor of the Bayes factor
	ln_pre_factor = gammaln(nu_c/2) + log (MAP_bb_prime * t_step * sqrt(pi) * sgn_bb_prime) - gammaln((1+nu_c)/2);
	
	% Integral factor of the Bayes factor
	z1 = - mu_c * sqrt(kappa_c / nu_c / sigma2_c);
	z2 = (MAP_bb_prime * t_step - mu_c) * sqrt(kappa_c / nu_c / sigma2_c);
	
	% Add 0 to integration points if it lies between the limits
	if z1 * z2 < 0
		z_int_points = [z1,0, z2];
	else
		z_int_points = [z1, z2];
	end;
	

	integrand_func = @(z) (1+z.^2).^(-(nu_c + 1)/2);

	% Split to increase accuracy
	hyperg = 0;
	for i = 1:length(z_int_points) - 1
		hyperg = hyperg + integral(integrand_func, z_int_points(i), z_int_points(i+1), 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE);
	end;
	ln_hyperg = log(hyperg * sgn_bb_prime);

	% Combine both into the Bayes factor K
	log_K(enum_conv_marginalized, bin) = (ln_pre_factor - ln_hyperg);
end;

log_K


figure(10);
clf;
hold on;

str_legend = {};
for convention = 1:conventions_count
	plot(x_bins_centers, log_K(convention, :), 'linewidth', line_width);
	str_legend{length(str_legend) + 1} = conventions_names{convention};
end;
x_lim_vec = xlim();

% Theory
plot(x_lim_vec, [0,0], 'k--', 'linewidth', line_width_theor);

ylim([-20,20]);
legend(str_legend);


















%% Calculate the theoretical bias in one bin, as expected based on the stationary solution of the Fokker-Planck equation
% The bias is calculated relative to supplied "true" values
%
% Parameters:
% bins_borders: (2, bins_count)
% lambda: scalar


function b_theor_bias = b_theor_bias_func(bins_borders, lambda, b_true)

%% Constants
load_constants;

%% Initialize
bins_count = size(bins_borders, 2);
% b_theor_bias = zeros(bins_count, 1);
b_mean = zeros(bins_count, 1);



%% Define the integrals used for <b> calculation
integrand_func = @(x, k) b_func(selected_D_case, x, L).^ k;
integral_func = @(x1, x2, k) integral(@(x) integrand_func(x, k), x1, x2, 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE);
b_mean_func = @(x1, x2, lambda) integral_func(x1, x2, 2 * lambda - 1) ./ integral_func(x1, x2, 2 * lambda - 2);



%% Calculate values for all bins
% <b> based on stationary Fokker-Planck
for bin = 1:bins_count
	b_mean(bin) = b_mean_func(bins_borders(1, bin), bins_borders(2, bin), lambda);
end;

% b bias
b_theor_bias = b_mean - b_true;





















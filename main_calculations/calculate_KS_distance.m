%% Calculate the Kolmogorov-Smirnov distance between a continuous and an empirical distribution
%
% This function does not use the algorithm proposed in Gonzalez et al., 1977,
% but a simpler algorithm also mentioned in that article


function distance = calculate_KS_distance(emp_sorted_data, cont_distrib_func)


%% Constants
ABS_TOLERANCE = 1e-3;
REL_TOLERANCE = 1e-3;

% % % % % Sort empirical distribution
% % % % emp_sorted_data = sort(emp_sorted_data);

% Initialize
n = length(emp_sorted_data);

% Define the continuous CDF
cum_cont_distrib_func = @(f) ... %integral(cont_distrib_func, -Inf, 0, 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE) +...
    integral(cont_distrib_func, -Inf, f, 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE);

% Evaluate CDF at empirical points
cum_cont_distrib = zeros(1, n);
for i = 1:n
    cum_cont_distrib(i) = cum_cont_distrib_func(emp_sorted_data(i));
	1;
end;



%% Calculate the distance with a simple algorithm
D_max_p = max((1:n) / n - cum_cont_distrib);
D_max_m = max(cum_cont_distrib - (0:n-1) / n);
distance = max(D_max_p, D_max_m);














%% Calculate the Kolmogorov-Smirnov distance between a continuous and an empirical distribution
%
% This function does not use the algorithm proposed in Gonzalez et al., 1977,
% but a simpler algorithm also mentioned in that article


function distance = calculate_KS_distance(emp_distrib, cont_distrib_func)


%% Constants
ABS_TOLERANCE = 1e-4;
REL_TOLERANCE = 1e-4;

% Sort empirical distribution
emp_distrib = sort(emp_distrib);
n = length(emp_distrib);

% Continuous CDF
cum_cont_distrib_func = @(f) integral(cont_distrib_func, -Inf, 0, 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE)...
    + integral(cont_distrib_func, 0, f, 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE);

% Evaluate CDF at empirical points
cum_cont_distrib = zeros(1, n);
for i = 1:n
    cum_cont_distrib(i) = cum_cont_distrib_func(emp_distrib(i));
end;



%% Calculate the distance with a simple algorithm
D_max_p = max((1:n) / n - cum_cont_distrib);
D_max_m = max(cum_cont_distrib - (0:n-1) / n);
distance = max(D_max_p, D_max_m);














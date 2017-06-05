


function out = calculate_KS_distance(emp_distrib, cont_distrib_func)


%% Constants
ABS_TOLERANCE = 1e-4;
REL_TOLERANCE = 1e-4;

% Sort empirical distribution
emp_distrib = sort(emp_distrib);
trials = length(emp_distrib);

% Continuous CDF
cum_cont_distrib_func = @(f) integral(cont_distrib_func, -Inf, 0, 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE)...
    + integral(cont_distrib_func, 0, f, 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE);

% Evaluate CDF at empirical points
cum_cont_distrib = zeros(1, trials);
for i = 1:trials
    cum_cont_distrib(i) = cum_cont_distrib_func(emp_distrib(i));
end;

%% Calculate distance with a simple algorithm
D_max_p = max((1:trials) / trials - cum_cont_distrib);
D_max_m = max(cum_cont_distrib - (0:trials-1) / trials);
out = max(D_max_p, D_max_m);














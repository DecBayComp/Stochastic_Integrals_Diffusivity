%% The function returns the location of the maximum and the distances to the
% 95% confidence intervals borders below and above



function out = find_confidence_interval(log_distr_func, max_search_range, bl_find_error_bars, MLE_guess, CONF_LEVEL, true_value)



%% Constants
ABS_TOLERANCE = 1e-7;
% ABS_LOG_TOLERANCE = 1e-5;
REL_TOLERANCE = 1e-7;
MAX_INITIAL_SEARCH_ITERATIONS = 100;
% max_range_to_search_for_root = 1e4;
bl_verbose = false;


%% Initialize
conf_percentile = (1 - CONF_LEVEL)/2;
distr_func = @(x) exp(log_distr_func(x));


%% Check if the function is 0 at the ends of the search interval. It should not be
% Reduce the search interval if it is
MLE_search_range = max_search_range;
half_width  = (max_search_range(2) - max_search_range(1)) / 2;
for i = 0:MAX_INITIAL_SEARCH_ITERATIONS
    if  log_distr_func(MLE_search_range(1)) > - Inf && log_distr_func(MLE_search_range(2)) > - Inf
        break;
    end;
    MLE_search_range = MLE_guess + [-1, 1] * half_width * (2/3) ^ i;
end;
if i > 0 
    fprintf('Search interval for MLE found after %i iterations.\n', i);
end;



%% Find the MLE value
neg_func_wrap = @(x) - log_distr_func(x);
optim_options = optimset('TolX', ABS_TOLERANCE);
MLE = fminbnd(neg_func_wrap, MLE_search_range(1), MLE_search_range(2), optim_options);
% If the function value at the MLE is not greater than 0, throw an error
if distr_func(MLE) <= 0
    error('The algorithm failed to find the MLE. Please modify the search region');
end;
if bl_verbose
    fprintf('MLE found: %.3f\n', MLE);
end;


%% Identify the scale of the function
% Find the drop by 'e' on one side




%% Find the confidence interval if requested
lower_boundary = MLE;
upper_boundary = MLE;
% Construct a cumulative function
cum_func = @ (x) integral(distr_func, -Inf, MLE, 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE)...
    + integral(distr_func, MLE, x, 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE);
if bl_find_error_bars
    tic;
    if bl_verbose
        fprintf('Calculating the confidence intervals for the current bin...\n');
    end;
    % Find a lower conf. interval boundary
    func_wrap = @(x) cum_func(x) - conf_percentile;
    optim_options = optimset('TolX', ABS_TOLERANCE, 'Display', 'notify');
    % Search for the initial interval
    interval = [max_search_range(1), MLE];
    interval = find_initial_interval_zero_search(func_wrap, interval, true);
    lower_boundary = ...
        fzero(func_wrap, interval, optim_options);
    if bl_verbose
        fprintf('Lower boundary found: %.3f\n', lower_boundary);
    end;
    % Find a higher conf. interval boundary
    func_wrap = @(x) cum_func(x) - (1 - conf_percentile);
    optim_options = optimset('TolX', ABS_TOLERANCE, 'Display', 'notify');
    % Search for the initial interval
    interval = [MLE, max_search_range(2)];
    interval = find_initial_interval_zero_search(func_wrap, interval, true);
    % Exact value
    upper_boundary = ...
        fzero(func_wrap, interval, optim_options);
    if bl_verbose
        fprintf('Upper boundary found: %.3f\n', upper_boundary);
        fprintf('FINISHED: Calculating the confidence intervals for the current bin. Time: %.2fs\n', toc);

        fprintf('Final result: %.3f <= %.3f <= %.3f\n', lower_boundary, MLE, upper_boundary);
    end;
    % disp([lower_boundary, MLE, higher_boundary]);
end;

%% Cacluate the minimum confidence interval that includes the true value
cum_prob = cum_func(true_value);
min_confidence_level = abs(1 - cum_prob * 2);



%% Prepare the output
out = [MLE, MLE - lower_boundary, upper_boundary - MLE, min_confidence_level];

1;






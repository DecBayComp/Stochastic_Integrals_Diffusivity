

function interval = find_initial_interval_zero_search(func, interval, bl_positive_derivative)


%% Constants
MAX_INITIAL_SEARCH_ITERATIONS = 1000;
MIN_INCREMENT = 1e-4;
C = 1e-2;


%% Initialize
%% Adjust function to take the derivative sign into account
if bl_positive_derivative
    func_growing = @(x) func(x);
else
    func_growing = @(x) - func(x);
end;
x_start = interval(1);
x_end = interval(2);
x_increment = max([abs(mean(interval)), MIN_INCREMENT]) * C;

% Left boundary
for i=0:MAX_INITIAL_SEARCH_ITERATIONS * 10
    if func_growing(x_start) > 0
        x_start = x_start - x_increment;
    else
        break;
    end;
end;
% Check convergence and throw error
if func_growing(x_start) > 0
    error('Unable to find the initial search interval for the error bars. This normally indicates that the max. search interval is too small. Consider extending');
end;
if i > 0 
    fprintf('Search interval found after %i iterations.\n', i);
end;
% Right boundary
for i=0:MAX_INITIAL_SEARCH_ITERATIONS
    if func_growing(x_end) < 0
        x_end = x_end + x_increment;
    else
        break;
    end;
end;
% Check convergence and throw error
if func_growing(x_end) < 0
    error('Unable to find the initial search interval for the error bars. This normally indicates that the max. search interval is too small. Consider extending');
end;
if i > 0 
    fprintf('Search interval found after %i iterations.\n', i);
end;

% Output interval
interval = [x_start, x_end];













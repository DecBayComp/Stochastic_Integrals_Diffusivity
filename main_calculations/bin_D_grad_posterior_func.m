
function result = bin_D_grad_posterior_func (l_ind, bin, D_j_grad, t_step)
%% Globals
% global SIGMA_MIN;
% global ABS_TOLERANCE;
% global REL_TOLERANCE;
global x_bins_centers_saved;
global D_MAX;

%% Constants
load_constants;

% 

%% Globals
% global x_bins_centers_saved;
% global V_j;

%% Initialization
if bin == 1
    disp('Error: cannot calculate gradient PDF for bin = 1. Skipping');
    result = zeros(size(D_j_grad));
    return;
end;
x_step = x_bins_centers_saved{l_ind}(bin) - x_bins_centers_saved{l_ind}(bin-1);
result = zeros(size(D_j_grad));
input_length = length(D_j_grad(:));


% Iterating over a gradient mesh
for index = 1:input_length
    cur_D_j_grad = D_j_grad(index);
    % Preparing the function to integrate
    integrand_func_wrap = @(Dj) x_step * bin_D_posterior_func (l_ind, bin-1, Dj - x_step * cur_D_j_grad, t_step) ...
        .* bin_D_posterior_func (l_ind, bin, Dj, t_step);
    
    % Determining the integration limits
    % {
    % Set the starting point of the integration
    D_limits = max([0, x_step * cur_D_j_grad]);
    

        %% Checking if peak points should be placed at the integration boundaries
        % The peak of the D distribtuion in bin (j-1)
        [~, ~, nu_n, sigma2_n] = get_n_parameters(l_ind, bin - 1);
        D_peak2 = nu_n * sigma2_n/(2*t_step * (2+nu_n)) + x_step * cur_D_j_grad;
        
        
        % The peak of the D distribtuion in bin j
        [~, ~, nu_n, sigma2_n] = get_n_parameters(l_ind, bin);
        D_peak1 = nu_n * sigma2_n/(2*t_step * (2+nu_n));
        
        
        if D_peak1 > D_limits(1)
            D_limits(2) = D_peak1;
        end;

        if D_peak2 > D_limits(1)
            D_limits(length(D_limits) + 1) = D_peak2;
        end;

        % Sorting the array of sigma points
        D_limits = sort(D_limits);

        % Determining the maximal sigma limit
        global_D_max = max([2 * max(D_limits), D_MAX]);
        D_limits(length(D_limits) + 1) = global_D_max;
    % }
    
    %% Integrating
    out = 0;
    for i = 1:(length(D_limits) - 1)
        out = out + integral(integrand_func_wrap, D_limits(i), D_limits(i+1), ...
            'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE);
    end;
    
    % Saving result for this particular gradient value
    result(index) = out;

end;

    






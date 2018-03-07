

function data_struct = infer_MAP_b(data_struct)

% Load constants
load_constants;

for bin = 1:x_bins_number
    %% Calculate MAP diffusivity D and b
    % Prepare functions
    log_D_function_to_minimze = @(D) bin_D_log_posterior_func (bin, D, t_step, data_struct, 'forward');
    log_b_function_to_minimze = @(b) bin_b_log_posterior_func (bin, b, t_step, data_struct, 'forward');
    
    % Make an MLE guess
    [mu_n, kappa_n, nu_n, sigma2_n] = get_n_parameters(bin, data_struct, 'forward');
    MLE_guess_D = 2 * nu_n / (nu_n + 2) * sigma2_n / (2 * t_step);
    MLE_guess_b = sqrt(MLE_guess_D);
    
    % Find confidence intervals if the bin is not empty
    if ~data_struct.bl_empty_bins(bin)
        D_inference = find_confidence_interval(log_D_function_to_minimze, [D_PRECISION, D_ABS_MAX], true, MLE_guess_D, CONF_LEVEL,...
            data_struct.D_theor_data(bin), trial, bin);
        b_inference = find_confidence_interval(log_b_function_to_minimze, [b_PRECISION, b_ABS_MAX], true, MLE_guess_b, CONF_LEVEL,...
            data_struct.b_theor_data(bin), trial, bin);
    else
        D_inference = ones(1, 4) * NaN;
        b_inference = ones(1, 4) * NaN;
    end
    
    % Save
    data_struct.MAP_D(bin, :) = D_inference;
    data_struct.MAP_b(bin, :) = b_inference;
end

 %% Regularize bb' gradient
% To obtain the bb' gradient, provide the function b^2/2 as input
b_squared_over_2 = data_struct.MAP_b(:, 1).^2 / 2;
[inferred_MAP_b_squared_over_2_reg, inferred_MAP_bb_prime_reg, inferred_MAP_bb_prime_reg_interpolated, norm_cost, x_grad_mesh] = ...
	regularize_gradient(b_squared_over_2, data_struct.x_bins_centers, alpha_reg);

% Save
data_struct.MAP_b_regular = sqrt(inferred_MAP_b_squared_over_2_reg * 2);
data_struct.MAP_bb_prime_regular = inferred_MAP_bb_prime_reg;
data_struct.MAP_bb_prime_regular_interp = inferred_MAP_bb_prime_reg_interpolated;
data_struct.x_grad_mesh = x_grad_mesh;


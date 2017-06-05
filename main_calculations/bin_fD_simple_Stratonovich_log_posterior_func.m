

function result = bin_fD_simple_Stratonovich_log_posterior_func(data_struct, l_ind, bin, fD, D_grad, str_direction)


%% Constants
load_constants;


%% Initialize
% Check grad length
if length(D_grad) > 1
    disp('Error: a scalar gradient variable required. Aborting');
    result = zeros(size(fD));
    return;
end;


%% Calculate
lambda = 1/2;
result = log(t_step / kBT) + bin_mu_log_posterior_func(data_struct, l_ind, bin, fD .* t_step ./ kBT + lambda * t_step * D_grad, str_direction);


















% % % % %% Globals
% % % % % global SIGMA_MIN;
% % % % % global ABS_TOLERANCE;
% % % % % global REL_TOLERANCE;
% % % % % global t_step;
% % % % global x_bins_centers_saved;
% % % % % global kBT;
% % % % 
% % % % %% Constants
% % % % load_constants;
% % % % % max_D = 10;
% % % % factor = 20;
% % % % 
% % % % %% Check bin number
% % % % % x_bins_number = length(x_bins_centers_saved{l_ind});
% % % % if bin == 1
% % % %     disp('Invalid bin number. Aborting');
% % % %     result = zeros(size(fD));
% % % %     return;
% % % % end;
% % % % 
% % % % 
% % % % %% Initialize
% % % % x_step = x_bins_centers_saved{l_ind}(bin) - x_bins_centers_saved{l_ind}(bin - 1);
% % % % 
% % % % 
% % % % %% Prepare D integration limits based on the maxima of posterior distributions (a good guess)
% % % % % bin j-1
% % % % [~, ~, nu_n1, sigma2_n1] = get_n_parameters(l_ind, bin - 1, data_struct);
% % % % % bin j
% % % % [~, ~, nu_n2, sigma2_n2] = get_n_parameters(l_ind, bin, data_struct);
% % % % 
% % % % % D1
% % % % D1_limits(1) = 0;
% % % % D1_limits(2) = nu_n1 * sigma2_n1 / (2 * t_step * (2 + nu_n1));
% % % % % Extend the region by the factor
% % % % D1_limits(3) = D1_limits(2) * factor;
% % % % 
% % % % 
% % % % % D2
% % % % D2_limits(1) = 0;
% % % % D2_limits(2) = nu_n2 * sigma2_n2 / (2 * t_step * (3 + nu_n2));
% % % % % Extend the region by the factor
% % % % D2_limits(3) = D2_limits(2) * factor;
% % % % 
% % % % 
% % % % 
% % % % %% Calculate for each value of fD independently
% % % % % fprintf('Starting integral calculation...\n');
% % % % input_length = length(fD);
% % % % result = zeros('like', fD);
% % % % for ind = 1:input_length
% % % %     tmp_result = 0;
% % % %     integrand_function_wrap = @(D1, D2) D1_D2_simple_Hanggi_integrand_function(l_ind, bin, D1, D2, fD(ind), t_step, x_step, kBT);
% % % %     for k1 = 1:2
% % % %         for k2=1:2
% % % %             tmp_result = tmp_result + integral2(integrand_function_wrap, D1_limits(k1), D1_limits(k1 + 1),...
% % % %                 D2_limits(k2), D2_limits(k2 + 1),...
% % % %                 'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE, 'Method', 'auto');
% % % %     %     result = result + integral2(integrand_function_wrap, 0, max_D, 0, @(D1) D1,...
% % % %     %         'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE, 'Method', 'auto');
% % % %     %     % toc;
% % % %     %     result = result + integral2(integrand_function_wrap, 0, max_D, @(D1) D1, max_D,...
% % % %     %         'RelTol', REL_TOLERANCE,'AbsTol', ABS_TOLERANCE, 'Method', 'auto');
% % % %         end;
% % % %     end;
% % % %     result(ind) = tmp_result;
% % % % end;
% % % % % fprintf('Calculation finished in %.1f s\n', toc);
% % % % 










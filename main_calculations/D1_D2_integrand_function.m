
function out = D1_D2_integrand_function(l_ind, bin, D1, D2, fD, t_step, x_step, kBT)


%% Initialize
[mu_n, kappa_n, nu_n, sigma2_n] = get_n_parameters(l_ind, bin);


%% Calculate
% 1. Error functions
x1 = sqrt(kappa_n ./ (4 .* D2 * t_step))...
    .* (fD * t_step / kBT - mu_n);
x2 = sqrt(kappa_n ./ (4 .* D2 .* t_step))...
    .* (fD * t_step / kBT + t_step .* (D2-D1) ./ x_step - mu_n);
%
log_out = log((erfc(x1) - erfc(x2)) ./ (D2-D1));

% zero_ind = find(out(:) == 0);
% log_out = log(out ./ (D2-D1));


% 1a. Adust for cases where D1==D2
ind = find(D1(:)==D2(:));
log_out(ind) = -(kappa_n * (fD * t_step - kBT * mu_n).^2 / kBT^2 + nu_n * sigma2_n)...
    ./ (4 * D2(ind) * t_step);


% 2). D-dependent coefficient
log_out = log_out - nu_n * sigma2_n ./ (4 * D2 * t_step)...
    + bin_sigma_squared_log_posterior_func (l_ind, bin - 1, 2 * D1 * t_step)...
    - (1 + nu_n/2) .* log(D2);

% 3. D-independent factor
log_out = log_out + (nu_n/2) * log(nu_n * sigma2_n / 4) - gammaln(nu_n/2)...
    + log(x_step) + (1 - nu_n/2) * log(t_step) - log(kBT);


% 4. Exponential
out = exp(log_out);


% % % % 5. Correct for zero difference of the error functions (to avoid infinity)
% % % out(zero_ind) = 0;

1;















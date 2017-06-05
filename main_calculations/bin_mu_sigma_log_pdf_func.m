

%% Defining the double distribution function for each individual bin
function out = bin_mu_sigma_log_pdf_func (l_ind, bin, mu, sigma2) 

%% Globals
% global n_j;
% global V_j;
% global dx_mean_in_bins_saved;

% % %% Reshaping mu and sigma arrays to make calculations on arrays
% % % Mu will be a column and sigma will be a row
% % mu_length = length(mu);
% % sigma_length = length(sigma);
% % 
% % % Reshaping in 1D
% % mu = reshape(mu, mu_length, 1);
% % sigma = reshape(sigma, 1, sigma_length);
% % 
% % % Converting to 2D
% % mu_2D = mu * ones(1, sigma_length);
% % sigma_2D = ones(mu_length, 1) * sigma;

%% Calculating
[mu_n, kappa_n, nu_n, sigma2_n] = get_n_parameters(l_ind, bin);
out = 1/2 * log(kappa_n / (2 * pi))...
        + nu_n / 2 * log(nu_n * sigma2_n / 2)...
        - gammaln(nu_n/2) - (nu_n + 3) / 2 .* log(sigma2)...
        - (kappa_n * (mu_n - mu).^2 + nu_n * sigma2_n) ./ (2 * sigma2);

% out = (...
%     (n_j{l_ind}(bin) + 1)/2 * log(n_j{l_ind}(bin)) - 1/2*log(pi)...
%     - (n_j{l_ind}(bin) - 1)/2 * log(2) - gammaln(n_j{l_ind}(bin)/2)...
%     + n_j{l_ind}(bin)/2 * log(V_j{l_ind}(bin))...
%     - (n_j{l_ind}(bin) + 2) .* log(sigma)...
%     - n_j{l_ind}(bin) / 2 ./ sigma.^2 ...
%     .* ((mu - dx_mean_in_bins_saved{l_ind}(bin)).^2 + V_j{l_ind}(bin))...
%     );


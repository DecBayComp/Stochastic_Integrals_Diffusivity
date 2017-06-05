

%% Defining the sigma distribution function for each individual bin
function out =  bin_sigma_log_likelihood_pdf_func (l_ind, bin, sigma)

%% Globals
global n_j;
global V_j;


alpha = (n_j{l_ind}(bin) -1) / 2;
beta = n_j{l_ind}(bin) * V_j{l_ind}(bin) / 2;

%% Calculating
out = log(n_j{l_ind}(bin)) .* (n_j{l_ind}(bin) - 1)./2 ...
    - log(2) .* (n_j{l_ind}(bin) - 3)./2 ...
    - gammaln((n_j{l_ind}(bin) - 1)./2)...
    + log(V_j{l_ind}(bin)) .* (n_j{l_ind}(bin) - 1)./2 ...
    - log(sigma) .* (n_j{l_ind}(bin)) ...
    - n_j{l_ind}(bin) .* V_j{l_ind}(bin) ./ 2 ./ sigma.^2;






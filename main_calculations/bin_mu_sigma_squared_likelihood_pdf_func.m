

function out = bin_mu_sigma_squared_likelihood_pdf_func (l_ind, bin, mu, sigma) 

out = exp(bin_mu_sigma_log_pdf_func(l_ind, bin, mu, sigma));
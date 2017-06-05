

function out = bin_sigma_likelihood_pdf_func (l_ind, bin, sigma) 

out = exp(bin_sigma_log_pdf_func(l_ind, bin, sigma));
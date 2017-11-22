%% Return the the part of data likelihood function to integrate over a, i.e. the a_part time the prior



function a_integrand = data_lklhd_cond_a_func_a_integrand(a, n, V, grad, mu_c, kappa_c, nu_c, sigma2_c, nu_pi, sigma2_pi, lambda, t_step, dx_mean_tot, V_tot)

[a_part, a_prior , ~, ~] = data_lklhd_cond_a_func(a, n, V, grad, mu_c, kappa_c, nu_c, sigma2_c, nu_pi, sigma2_pi, lambda, t_step, dx_mean_tot, V_tot);
a_integrand = a_part .* a_prior;























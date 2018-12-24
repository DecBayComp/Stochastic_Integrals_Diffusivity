
try:
    has_run
except NameError:
    %matplotlib
    %load_ext autoreload
    %autoreload 2

    has_run = 1
else:
    print("Graphic interface NOT re-initialized")

# import this

import numpy as np
from scipy.integrate import dblquad, quad
from scipy.special import factorial, gammainc, gammaln

from tramway.inference.bayes_factors.calculate_bayes_factors import \
    calculate_bayes_factors
from tramway.inference.bayes_factors.calculate_marginalized_integral import (calculate_integral_ratio,
                                                                             calculate_marginalized_integral)
from tramway.inference.bayes_factors.calculate_posteriors import (calculate_one_1D_posterior_in_2D,
                                                                  calculate_one_1D_prior_in_2D,
                                                                  calculate_one_2D_posterior)
from tramway.inference.bayes_factors.convenience_functions import n_pi_func, p

# %%
zeta_t = np.array([-0.008164965809277261, 0.0])
zeta_sp = np.array([-0.01632993, -0.])
n = 20
u = 0.95
dim = 2
n_pi = n_pi_func(dim)
# rel_loc_error = 1.47
loc_error = 0.1**2.0
eta = np.sqrt(n_pi / (n + n_pi))
pow = p(n, dim)
v0 = 1 + n_pi / n * u
V = 0.0024
V_pi = V * u

# res = calculate_marginalized_integral(zeta_t=zeta_t, zeta_sp=zeta_sp,
#                                       p=pow, v=v0, E=eta ** 2.0, rel_loc_error=rel_loc_error, zeta_a = zeta_a)


# calculate_one_posterior(zeta_t = zeta_t, zeta_sp=zeta_sp, zeta_a=zeta_a, n = n, V=V, V_pi=V_pi, loc_error= loc_error, dim=2, lamb='marg')
#
#
# # Check norm
# def integrate_me(zay, zax):
#     return calculate_one_posterior(zeta_t = zeta_t, zeta_sp=zeta_sp, zeta_a=[zax, zay], n = n, V=V, V_pi=V_pi, loc_error= loc_error, dim=2, lamb='marg')
#
# lims = np.array([-1,1])*5
# tol = 1e-4
# r1= dblquad(integrate_me, lims[0], lims[1], lims[0], lims[1], epsabs=tol, epsrel=tol)
# print(r1)

zeta_sp = [0.1, 0]
zeta_t = [0, 0]
zeta_a = 0
V = np.asarray(0.8 ** 2.0)
u = np.asarray(0.95)
V_pi = u * V
loc_error = 0.1**2.0
n = 500
calculate_one_1D_prior_in_2D(zeta_a=zeta_a, V_pi=V_pi, loc_error=loc_error)


# calculate_one_1D_posterior_in_2D(zeta_t=zeta_t, zeta_sp=zeta_sp, zeta_a=zeta_a,
#                                  n=n, V=V, V_pi=V_pi, loc_error=loc_error, lamb='marg', axis='x')

# # %% Check norm
#
#
# def integrate_me(zax):
#     return calculate_one_1D_posterior_in_2D(zeta_t=zeta_t, zeta_sp=zeta_sp, zeta_a=zax, n=n, V=V, V_pi=V_pi, loc_error=loc_error, lamb=0.1)
#
#
# lims = np.array([-1, 1]) * 10
# tol = 1e-4
# r1 = quad(integrate_me, lims[0], lims[1], epsabs=tol, epsrel=tol)
# print(r1)
2.**400
1.007125**(-401)
factorial(401)

gammainc(1, 0)

exp(1)
np.log(0)

np.inf == True
1 == True

if np.inf:
    print(1)

1 / np.array(0)

zeta_t = np.array([1, 1])
zeta_sp = np.array([1, 0])
with np.errstate(invalid='ignore', divide='ignore'):
    x = zeta_t / zeta_sp
x

a, b = 1, 2

1 / 0.96


# %%
zeta_t = np.asarray([0.7, 0.4])
zeta_sp = np.asarray([0.8, 0.6])
n = np.asarray(20)
V = np.asarray(0.8 ** 2.0)
u = np.asarray(0.95)
loc_error = 0.1**2.0
V_pi = u * V
lamb = 0  # 'marg'
tol = 1e-4
zax = 0

# Wrapper to get just a function of zeta_ax, zeta_ay
calculate_one_1D_posterior_in_2D(zeta_t=zeta_t, zeta_sp=zeta_sp,
                                 zeta_a=zax, n=n, V=V, V_pi=V_pi, loc_error=loc_error, lamb=lamb)


# >> Test zero localization error <<
zeta_ts = np.asarray([[0.7, 0.4]])
zeta_sps = np.asarray([[0.8, 0.6]])
ns = np.asarray([[20]])
Vs = np.asarray([[0.8 ** 2.0]])
us = np.asarray([[0.95]])
loc_error = 0.1**2.0 * 0
Vs_pi = us * Vs
lg_Bs, forces, _ = calculate_bayes_factors(
    zeta_ts=zeta_ts, zeta_sps=zeta_sps, ns=ns, Vs=Vs, Vs_pi=Vs_pi, loc_error=loc_error)
lg_Bs
np.exp(lg_Bs)


# %% >> Test Bayes factor in one input bin <<
zeta_ts = np.asarray([[0.7, 0.4]])
zeta_sps = np.asarray([[0.8, 0.6]])
ns = np.asarray([[20]])
Vs = np.asarray([[0.8 ** 2.0]])
us = np.asarray([[0.95]])
loc_error = 0
Vs_pi = us * Vs
lg_Bs, forces, _ = calculate_bayes_factors(
    zeta_ts=zeta_ts, zeta_sps=zeta_sps, ns=ns, Vs=Vs, Vs_pi=Vs_pi, loc_error=loc_error)
print(7, lg_Bs, 10**lg_Bs)
true_B = 0.2974282533

# # Check value
# self.assertTrue(np.isclose(10**lg_Bs[0], true_B, rtol=self.rel_tol, atol=self.tol),
#                 "Bayes factor calculation failed for one bin. The obtained B = %.8g does not match the expected B = %.8g" % (10**lg_Bs[0, 0], true_B))
# # Check force presence
# self.assertTrue((true_B >= self.B_threshold) ==
#                 forces[0], "Boolean conservative force return incorrect for the case of one bin")
# 0.09 + 0.01

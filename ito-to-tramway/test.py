
try:
    has_run
except NameError:
    %matplotlib
    %load_ext autoreload
    %autoreload 2

    has_run = 1
else:
    print("Graphic interface NOT re-initialized")

import matplotlib.pyplot as plt
import numpy as np

from tramway.inference.bayes_factors.calculate_bayes_factors import (_calculate_one_bayes_factor,
                                                                     calculate_minimal_n)
from tramway.inference.bayes_factors.get_D_posterior import (get_D_confidence_interval,
                                                             get_D_posterior,
                                                             get_MAP_D)

# import scipy.integrate as integrate


# %
n = 121
zeta_t = [-0.11129124,  0.12443088]
V = 0.002323199914280775
V_pi = 0.0023187586103097994
sigma2 = 1e-8
dim = 2
dt = 0.04

D_step = 1e-4
D_max = 0.03
D_mesh = np.arange(0, D_max + D_step, D_step)


get_MAP_D(n=n, zeta_t=zeta_t, V=V, V_pi=V_pi, dt=dt, sigma2=sigma2, dim=dim)

post_func = get_D_posterior(n=n, zeta_t=zeta_t, V=V, V_pi=V_pi, dt=dt, sigma2=sigma2, dim=dim)
post_func(0.297)


# %%
post_data = [post_func(D) for D in D_mesh]
# check_norm = integrate.quad(post_func, 0, 1e6)
norm = np.trapz(post_data, D_mesh)
print(norm)

# post_data


# %% Confidence intervals
alpha = 0.95
MAP_D, CI = get_D_confidence_interval(alpha=alpha, n=n, zeta_t=zeta_t, V=V,
                                      V_pi=V_pi, dt=dt, sigma2=sigma2, dim=dim)
print(MAP_D, CI)


# %% Plot
color = 'rebeccapurple'
fig, ax = plt.subplots(num=1, clear=True)
ax.plot(D_mesh, post_data)
ylims = plt.ylim()

for x in CI:
    ax.plot([x] * 2, ylims, '--', color=color)

ax.plot([MAP_D] * 2, ylims, 'k--')


fig.show()

# %% Test minimal n calculation
calculate_minimal_n(zeta_t=np.array([np.nan, np.nan]), zeta_sp=np.array([np.nan, np.nan]), n0=1,
                    V=np.nan, V_pi=1, loc_error=1e-8, dim=2, B_threshold=10)
_calculate_one_bayes_factor(zeta_t=None, zeta_sp=np.array([np.nan, np.nan]), n=11,
                            V=np.nan, V_pi=1, loc_error=1e-8, dim=2, B_threshold=10)


np.nan is None


# import numpy as np
# from scipy.integrate import dblquad, quad
# from scipy.special import factorial, gammainc, gammaln
#
# from calculate import calculate
# from constants import dt
# # import this
# from tesselate_and_infer import tesselate_and_infer
#
# # from tramway.inference.bayes_factors.calculate_bayes_factors import \
# #     calculate_bayes_factors
# # from tramway.inference.bayes_factors.calculate_marginalized_integral import (calculate_integral_ratio,
# #                                                                              calculate_marginalized_integral)
# # from tramway.inference.bayes_factors.calculate_posteriors import (calculate_one_1D_posterior_in_2D,
# #                                                                   calculate_one_1D_prior_in_2D,
# #                                                                   calculate_one_2D_posterior)
# # from tramway.inference.bayes_factors.convenience_functions import n_pi_func, p
#
# # %% Produce a diffusivity map for one of the simulated trajectories for the article
# file = r'D:\Google Drive\git\Stochastic_Integrals_Diffusivity\ito-to-tramway\input\diffusivity_map_for_article_sim_trajectory\sim_data_000000050.csv'
# tesselate_and_infer(file, localization_error=0, load=False)
#
# output_folder = r'D:\Google Drive\git\Stochastic_Integrals_Diffusivity\ito-to-tramway\input\diffusivity_map_for_article_sim_trajectory\result'
# calculate(file, output_folder, bl_produce_maps=True,
#           snr_label='snr', localization_error=0, ticks=True)

# %% Magic functions
try:
    has_run
except NameError:
    %matplotlib
    %load_ext autoreload
    %autoreload 2

    has_run = True
#

import importlib
import sys

import numpy as np
from tqdm import tqdm, trange

from calculate_bayes_factor_limits import calculate_bayes_factor_limits
# %% imports
from calculate_zeta_t_lims import calculate_zeta_t_lims
from calculate_zeta_t_roots import calculate_zeta_t_roots
from find_marginalized_zeta_t_roots import find_marginalized_zeta_t_roots
from log_C import log_C as log_C_func
from plot_K_L import plot_K_L
from plot_K_L_2D import plot_K_L_2D
from plot_zeta_t_perp import plot_zeta_t_perp

# %% Main calculations

# Constants
lambs = (0, 0.5, 1.0)
zeta_sp_abs_lim = 1.00  # 0.02
zeta_sp_steps = 100
dim = 2
n_pi = 5 - dim  # minimum number of jumps is different for priors in 1D and 2D
B = 1  # level of evidence
log10_B = np.log10(B)

lambs_count = len(lambs)

# Define dimension-specific parameters
# Relative prior uncertainties (u = Vp/V)
us_1D = [0.1, 10.0]
us_2D = [0.1, 10.0]

# total force component orthogonal to the diffusivity gradient
zeta_t_perp_1D = [0.0]
zeta_t_perp_2D = [0.0]

# number of jumps
ns_1D = [5, 20, 100]
ns_2D = [5, 20, 100]

#
if dim == 1:
    us = us_1D
    ztpers = zeta_t_perp_1D
    ns = ns_1D
else:
    us = us_2D
    ztpers = zeta_t_perp_2D
    ns = ns_2D

us_count = len(us)
ztpers_count = len(ztpers)
ns_count = len(ns)
zeta_t_roots_size = (us_count, ns_count, lambs_count + 1, 2, zeta_sp_steps)

# zeta_t_roots_size = (ztpers_count, ns_count,
#                      lambs_count + 1, 2, zeta_sp_steps)

# Create zeta_sp mesh
zeta_sps = np.linspace(- zeta_sp_abs_lim, zeta_sp_abs_lim, zeta_sp_steps)

zeta_t_roots = np.zeros(zeta_t_roots_size, dtype=np.float) * np.nan


#     v = 1 + u * n_pi / ns + (dim - 1) * s * zeta_t_perp ** 2

for n_ind in range(ns_count):
    n = ns[n_ind]
    # print("Processing u: %i/%i, n: %i/%i" %
    #       (u_ind + 1, us_count, n_ind + 1, ns_count))

    for u_ind in range(us_count):
        u = us[u_ind]

        for ztper_ind in range(ztpers_count):
            ztper = ztpers[ztper_ind]

            # Check if the requested Bayes factor can be reached
            log10_B_lims = calculate_bayes_factor_limits(
                ns=ns, dim=dim, zeta_t_perp=ztper, u=u)
            bl_reachable = log10_B >= log10_B_lims[n_ind,
                                                   0] and log10_B <= log10_B_lims[n_ind, 1]

            if not bl_reachable:
                print("Warning: The value of the Bayes factor (B = %.3g) cannot be achieved for the current parameters. Bayes factor limits: [%.3g; %.3g]"
                      % (10 ** log10_B, 10 ** log10_B_lims[n_ind, 0], 10 ** log10_B_lims[n_ind, 1]))
                continue

            for zeta_ind in trange(zeta_sp_steps):
                zeta_sp = zeta_sps[zeta_ind]

                # Fixed-lambda
                for lamb_ind in range(lambs_count):
                    lamb = lambs[lamb_ind]

                    cur_roots = calculate_zeta_t_roots(
                        zeta_sp=zeta_sp, lamb=lamb, n=n, B=B, dim=dim, zeta_t_perp=ztper, u=u)

                    # Save
                    # if dim == 1:
                    zeta_t_roots[u_ind, n_ind, lamb_ind,
                                 0, zeta_ind] = cur_roots[0]
                    zeta_t_roots[u_ind, n_ind, lamb_ind,
                                 1, zeta_ind] = cur_roots[1]
                    # else:
                    #     zeta_t_roots[ztper_ind, n_ind, lamb_ind,
                    #                  0, zeta_ind] = cur_roots[0]
                    #     zeta_t_roots[ztper_ind, n_ind, lamb_ind,
                    #                  1, zeta_ind] = cur_roots[1]

                # Marginalized
                cur_roots = find_marginalized_zeta_t_roots(
                    zeta_sp=[zeta_sp], n=n, n_pi=n_pi, B=B, u=u, dim=dim, zeta_t_perp=ztper)

                # Save
                # if dim == 1:
                zeta_t_roots[u_ind, n_ind, lambs_count,
                             0, zeta_ind] = cur_roots[0]
                zeta_t_roots[u_ind, n_ind, lambs_count,
                             1, zeta_ind] = cur_roots[1]
                # else:
                #     zeta_t_roots[ztper_ind, n_ind, lambs_count,
                #                  0, zeta_ind] = cur_roots[0]
                #     zeta_t_roots[ztper_ind, n_ind, lambs_count,
                #                  1, zeta_ind] = cur_roots[1]

print("Calculation terminated")
# print(zeta_t_roots)

# print(zeta_t_roots)

# %% plot
# if dim == 1:
plot_K_L(zeta_sps=zeta_sps, zeta_t_roots=zeta_t_roots,
         ns=ns, us=us, dim=dim, zeta_t_perp=ztpers)
# else:
#     plot_K_L(zeta_sps=zeta_sps, zeta_t_roots=zeta_t_roots,
#              ns=ns, us=us, dim=dim, zeta_t_perp=ztpers)


#

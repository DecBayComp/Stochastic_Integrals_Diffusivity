"""
Estimate expected statistical performance of the active force detection criterion
"""

from constants import D_0, D_ratio, k, dt
from bayes_factors.find_marginalized_zeta_t_roots import find_marginalized_zeta_t_roots
import numpy as np


def estimate_theoretical_performance(data, expect_mean_n):
    dim = 2
    n_pi = 5 - dim
    u = 1.0
    sim_zeta_t_y_over_zeta_sp = 6.25
    theor_Bs = [0.1, 10.0]
    D_max = D_0 * D_ratio
    sim_zeta_sp_par = k * D_0 * np.sqrt(dt) / (np.sqrt(D_0) + np.sqrt(D_max))
    print(sim_zeta_sp_par)

    # Zeta_t roots calculation (probably needs some testing)
    exp_zeta_ts = np.zeros((2, 4), dtype=np.float32)
    exp_zeta_ts_over_zeta_sps = np.zeros((2, 4), dtype=np.float32)

    # Low B threshold for no perp
    zeta_t_perps = np.asarray([0.0, sim_zeta_t_y_over_zeta_sp]) * sim_zeta_sp_par
    print("Expected zt_per/zsp: %.3f" % (sim_zeta_t_y_over_zeta_sp * sim_zeta_sp_par))
    for f_ind in range(2):
        zeta_t_perp = zeta_t_perps[f_ind]
        # print(zeta_t_perp)

        exp_zeta_ts[f_ind, [1, 2]] = find_marginalized_zeta_t_roots(zeta_sp_par=sim_zeta_sp_par,
                                                                    n=expect_mean_n[f_ind], n_pi=n_pi, B=theor_Bs[0], u=u, dim=dim, zeta_t_perp=zeta_t_perp)
        exp_zeta_ts[f_ind, [0, 3]] = find_marginalized_zeta_t_roots(zeta_sp_par=sim_zeta_sp_par,
                                                                    n=expect_mean_n[f_ind], n_pi=n_pi, B=theor_Bs[1], u=u, dim=dim, zeta_t_perp=zeta_t_perp)

    exp_zeta_ts_over_zeta_sps = exp_zeta_ts / sim_zeta_sp_par
    print("Finished!\nzt_par / zsp roots found:")
    print(exp_zeta_ts_over_zeta_sps)

    return exp_zeta_ts_over_zeta_sps

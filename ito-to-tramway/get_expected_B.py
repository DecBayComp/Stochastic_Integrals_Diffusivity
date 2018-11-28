"""
The function returns a function which calculates the theoretical Bayes factor given ksi and ztpar as input. The other experimental parameters are implicitly loaded here.

    Comments:
    It doesn't matter if zsp has a positive projection on x, as long as I keep the ratio ztx / zsp of the right sign.

    # TODO: check that the variance is correctly estimated. I think it misses sqrt(2)
"""

import numpy as np

from constants import D_0, D_ratio, L, dt, zeta_t_y_over_zeta_sp_abs
from tramway.inference.bayes_factors import calculate_bayes_factors


def get_expected_B(ksis, D_sims, ns, bl_ztpers):
    # Change these
    u = 1.0
    dim = 2

    # Do not change these
    V = 1
    V_pi = V * u
    D_grad_abs = 2 * D_0 / L * (D_ratio - 1)  # um^2/s
    # print(D_grad)
    M = len(ksis)
    lg_Bs = np.ones(M) * np.nan

    # Convert to format accepted by the calculation function
    zeta_sps = np.zeros((M, 2))
    zeta_ts = np.zeros((M, 2))
    Vs = np.ones(M)
    zeta_sps[:, 0] = D_grad_abs * np.sqrt(dt / 2 / D_sims / dim)

    zeta_ts[:, 0] = zeta_sps[:, 0] * ksis
    zeta_ts[:, 1] = zeta_sps[:, 0] * zeta_t_y_over_zeta_sp_abs * bl_ztpers
    Vs_pi = Vs * u

    # Calculate the Bayes factor
    lg_Bs, _, _ = calculate_bayes_factors(zeta_ts=zeta_ts, zeta_sps=zeta_sps, ns=ns,
                                          Vs=Vs, Vs_pi=Vs_pi, loc_error=0, dim=2, verbose=True)

    return lg_Bs

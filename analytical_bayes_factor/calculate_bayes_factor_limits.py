

import numpy as np


def calculate_bayes_factor_limits(ns, n_pi, dim, zeta_t_perp, u):
    """
    The function calculates limits of achievable Bayes factor and outputs them as base 10 logarithms.
    """

    ns = np.transpose(np.asarray(ns))
    etas = np.sqrt(n_pi / (ns + n_pi))
    ps = dim * (ns + n_pi - 1.0) / 2.0 - 1.0

    # Define v function
    def vs(s):
        return(1 + n_pi / ns * u + (dim - 1) * s * zeta_t_perp ** 2)

    # print(vs(1))
    low_lims = np.multiply(ps, np.log10(vs(1))) - \
        np.multiply(ps, np.log10(vs(etas ** 2)))
    low_lims = low_lims + dim * np.log10(etas)

    upper_lims = np.multiply(-2 * ps + dim, np.log10(etas))

    log10_B_lims = np.stack((low_lims, upper_lims), axis=1)
    # print(log10_B_lims)
    # print(10**log10_B_lims)

    return log10_B_lims

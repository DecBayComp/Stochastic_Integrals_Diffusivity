

import numpy as np
from scipy import integrate


def calculate_marginalized_integral(zeta_t, zeta_sp, p, v, E, dim):
    """
    Calculate the lambda integral

    Integrate[(v + E * (zeta_t - lambda * zeta_sp) **2) ** (-p), {lambda, 0, 1}]

    for the given values of zeta_t and zeta_sp.
    They must be 1D lists of same length (supposedly, treated as x projections for different parameter values).
    n and n0 are scalars.
    """

    # Constants
    tol = 1e-8

    # Check input
    if not isinstance(zeta_t, list) or not isinstance(zeta_sp, list):
        raise TypeError("'zeta_t' and 'zeta_sp' must be 1D lists.")

    zeta_length = len(zeta_t)
    if zeta_length != len(zeta_sp):
        raise ValueError("'zeta_t' and 'zeta_sp' must have same length.")

    results = []
    for zeta_ind in range(zeta_length):
        zeta_sp_cur = zeta_sp[zeta_ind]
        zeta_t_cur = zeta_t[zeta_ind]

        # If zeta_sp ~ 0, no need to actually integrate
        if np.isclose(zeta_sp_cur, 0, rtol=tol, atol=tol):
            result = (v + E * zeta_t_cur ** 2.0) ** (-p)
            results.append([result, 0])
            continue

        # Check if a break-point needs to be added
        break_lambda = zeta_t_cur / zeta_sp_cur
        if break_lambda < 0.0 or break_lambda > 1.0:
            break_lambda = None
        else:
            break_lambda = [break_lambda]

        # Define the integrand function
        def integrate_me(l):
            return (v + E * (l * zeta_sp_cur - zeta_t_cur) ** 2.0) ** (-p)

        # Perform integration
        result = integrate.quad(integrate_me, 0.0, 1.0, points=break_lambda, full_output=0)

        results.append(result)

    return(results)

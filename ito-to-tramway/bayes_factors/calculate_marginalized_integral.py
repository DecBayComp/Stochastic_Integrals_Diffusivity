# Copyright © 2018, Alexander Serov

import numpy as np
from scipy import integrate
from scipy.special import gamma, gammainc, gammaincc


def calculate_marginalized_integral(zeta_t, zeta_sp, p, v, E, rel_loc_error, abs_tol=1e-8):
    """
    Calculate the marginalized lambda integral
    >>>
    Integrate[gamma_inc[p, arg * rel_loc_error] * arg ** (-p), {lambda, 0, 1}]
    >>>
    for the given values of zeta_t and zeta_sp.

    Here:

    arg = (v + E * (zeta_t - lambda * zeta_sp)**2);
    gamma_inc --- is the non-normalized lower incomplete gamma function;
    rel_loc_error = n * V / (2 * dim * sigma_L^2) --- inverse relative localization error.

    Input:
    zeta_t and zeta_sp are vectors of length 2 (x, y).
    All other parameters are scalars.

    abs_tol --- absolute tolerance for integral calculations.
    """

    # Constants
    abs_tol = 1e-8
    rel_tol = 1e-8

    def gammainc_(a, x):
        """Non-normalized incomplete gamma function."""
        return gammainc(a, x)  # * gamma(a)

    # Convert to numpy
    zeta_t = np.asarray(zeta_t)
    zeta_sp = np.asarray(zeta_sp)

    def diff(l):
        return l * zeta_sp - zeta_t

    def arg(l):
        return v + E * (diff(l) @ diff(l).T)

    # If norm(zeta_sp) ~ 0, no need to actually integrate
    if np.linalg.norm(zeta_sp) < abs_tol:
        result = gammainc_(p, arg(0) * rel_loc_error) * arg(0) ** (-p)
        return (result)

    # Check if break points need to be added
    # ignore if 0 is divided by 0
    with np.errstate(invalid='ignore', divide='ignore'):
        lambda_breaks_candidates = np.divide(zeta_t, zeta_sp)
    lambda_breaks = [lb for lb in lambda_breaks_candidates if lb >= 0.0 and lb <= 1.0]

    # Define the integrand function for cases when the localization error is negligible and otherwise
    def get_integrate_me():

        def with_gamma(l):
            # print((p, arg(l) * rel_loc_error, gammainc(p, arg(l) * rel_loc_error),
            #        gammaincc(p, arg(l) * rel_loc_error)))
            return gammainc_(p, arg(l) * rel_loc_error) * arg(l) ** (-p)

        def without_gamma(l):
            return arg(l) ** (-p)

        assert rel_loc_error >= 0, f"Negative relative localization error: {rel_loc_error}."
        if 1 / rel_loc_error <= abs_tol:
            return without_gamma
        else:
            return with_gamma

    integrate_me = get_integrate_me()

    # Perform integration
    result = integrate.quad(integrate_me, 0.0, 1.0,
                            points=lambda_breaks, full_output=0, epsabs=abs_tol, epsrel=rel_tol)
    # print(result)

    return(result[0])

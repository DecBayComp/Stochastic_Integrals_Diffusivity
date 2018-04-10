


from log_C import log_C
import numpy as np


def ito_roots_equation(zeta_t, n_0, n, n_pi):
	"""
	Calculate function vaue that defines zeta_t roots for the lambda = 0 case.
	Zeta_t may be a vector
	"""
	zeta_t = np.asarray(zeta_t)

	res = - 2.0 / (n_0 - 2.0) * log_C(n_0, n_pi)
	res = np.exp(res)
	par = (1 + n * n_pi / n_0 ** 2.0 * zeta_t ** 2.0) ** ((n_0 - 3.0) / (n_0 - 2.0))
	res = res * par - 1.0 - n / n_0 * zeta_t ** 2.0
	return(res)


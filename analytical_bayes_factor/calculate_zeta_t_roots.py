

from log_C import log_C
import numpy as np
from scipy.special import gammaln	# for gammaln

def calculate_zeta_t_roots(zeta_sp, lamb, n, n_pi, K):
	"""
	Calculate zeta_t roots for the given parameters.
	All inputs should be scalar.
	Output: [zeta_t-, zeta_t+]
	"""
	n0 = n + n_pi

	zeta_t_roots = 2.0 / (n0 - 2.0) * (np.log(K) - log_C(n0, n_pi))
	zeta_t_roots = np.exp(zeta_t_roots) - 1.0
	zeta_t_roots = np.sqrt(n / n0 * zeta_t_roots)
	zeta_t_roots = lamb * zeta_sp + np.asarray([-1.0, 1.0]) * zeta_t_roots

	return zeta_t_roots

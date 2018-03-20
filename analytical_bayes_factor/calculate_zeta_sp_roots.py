

from log_C import log_C
import numpy as np
from scipy.special import gammaln	# for gammaln

def calculate_zeta_sp_roots(zeta_t, lamb, n_0, n, n_pi):
	"""
	Calculate zeta_sp roots for the given parameters and K_L = 1.
	All inputs should be scalar.
	Output: [zeta_sp+, zeta_sp_-]
	"""

	zeta_sp_roots = -2.0 / (n_0 - 2.0) * log_C(n_0, n_pi)
	par = (1 + n * n_pi / n_0 ** 2.0 * zeta_t ** 2.0) ** ((n_0 - 3.0) / (n_0 - 2.0))
	zeta_sp_roots = np.exp(zeta_sp_roots) * par - 1.0
	# print(zeta_sp_roots)
	zeta_sp_roots = np.sqrt(zeta_sp_roots * n_0 / n)
	zeta_sp_roots = (zeta_sp_roots * np.asarray([-1.0, 1.0]) + zeta_t) / lamb

	return zeta_sp_roots

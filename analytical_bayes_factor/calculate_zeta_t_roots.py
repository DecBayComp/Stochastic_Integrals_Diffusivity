

from log_C import log_C
import numpy as np
from scipy.special import gammaln	# for gammaln

def calculate_zeta_t_roots(zeta_sp, lamb, n, n_pi, B, v):
	"""
	Calculate zeta_t roots for the given parameters.
	All inputs should be scalar.
	Output: [zeta_t-, zeta_t+]
	"""
	eta = n_pi / (n + n_pi)

	upstairs = (np.log10(B) - np.log10(eta)) * 2.0 / (n + n_pi - 3)
	upstairs = 1.0 - 10.0 ** upstairs

	downstairs = (np.log10(B) - np.log10(eta)) * 2.0 / (n + n_pi - 3)
	downstairs += 2.0 * np.log10(eta)
	downstairs = 10.0 ** downstairs - 1

	zeta_t_roots = v * upstairs / downstairs
	zeta_t_roots = np.sqrt(zeta_t_roots)
	zeta_t_roots = lamb * zeta_sp + np.asarray([-1.0, 1.0]) * zeta_t_roots
	# print(zeta_t_roots)

	return zeta_t_roots

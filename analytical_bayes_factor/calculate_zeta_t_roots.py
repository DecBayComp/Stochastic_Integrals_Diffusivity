

from log_C import log_C
import numpy as np
from scipy.special import gammaln	# for gammaln

def calculate_zeta_t_roots(zeta_sp, lamb, n, n_pi, B, dim, zeta_t_perp, u):
	"""
	Calculate zeta_t roots for the given parameters.
	All inputs should be scalar.
	Output: [zeta_t-, zeta_t+]
	"""

	eta = np.sqrt(n_pi / (n + n_pi))
	p = dim * (n + n_pi + 1.0) / 2.0 - 2.0

	# Define v function
	def v(s):
		return(1.0 + n_pi / n * u + (dim - 1.0) * s * zeta_t_perp ** 2)

	upstairs = (np.log(B) - dim * np.log(eta)) / p
	upstairs = v(1.0) - v(eta ** 2.0) * np.exp(upstairs)

	downstairs = (np.log(B) - dim * np.log(eta)) / p
	downstairs += 2.0 * np.log(eta)
	downstairs = np.exp(downstairs) - 1.0

	zeta_t_roots = upstairs / downstairs
	zeta_t_roots = np.sqrt(zeta_t_roots)
	zeta_t_roots = lamb * zeta_sp + np.asarray([-1.0, 1.0]) * zeta_t_roots
	# print(zeta_t_roots)

	return zeta_t_roots

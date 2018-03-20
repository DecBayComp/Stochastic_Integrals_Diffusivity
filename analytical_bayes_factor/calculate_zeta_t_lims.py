

import numpy as np
from log_C import log_C
from scipy.special import gammaln	# for gammaln


def calculate_zeta_t_lims(n_0, n, n_pi):
	"""
	Returns [zeta_t_lims, bl_lims_exist], where `bl_lims_exist` indicates whether the limits exist. 
	Only then the return values in `zeta_t_lims` are meaningful. 
	Otherwise the requested Bayes factor exists in the whole complex plain, with no limits
	"""

	zeta_t_lims = 2.0 / (n_0 - 3.0) * log_C(n_0, n_pi)
	zeta_t_lims = np.exp(zeta_t_lims) - 1.0

	# If zeta_t_lims is now negative, there are no borders.
	# Physically this means that this Bayes factor value can be achieved everywhere in the region
	if zeta_t_lims < 0:
		print("Zeta_t limits do not exist for these points number and K_L")
		return [[float('nan'), float('nan')], False]

	zeta_t_lims = np.sqrt(n_0 ** 2.0 / n / n_pi * zeta_t_lims) * np.asarray([-1.0, 1.0])
	# print (zeta_t_lims)
	return [zeta_t_lims, True]


import numpy as np
from scipy.special import gammaln	# for gammaln


def calculate_zeta_t_lims(n_0, n, n_pi):
	"""
	Returns [zeta_t_lims, bl_lims_exist], where `bl_lims_exist` indicates whether the limits exist. Only then the return values in `zeta_t_lims` are meaningful. Otherwise the requested Bayes factor does not exist for this `zeta_t` and `n`
	"""

	zeta_t_lims =  gammaln((n_pi - 2) / 2) + gammaln((n_0 - 3)/2) - gammaln((n_pi - 3) / 2) - gammaln((n_0 - 2) / 2) + (n_pi - 2) / 2 * np.log(n_0 / n_pi)
	zeta_t_lims = 2 / (n_0 - 3) * zeta_t_lims
	zeta_t_lims = np.exp(zeta_t_lims) - 1.0

	# If zeta_t_lims is now negative, there are no borders.
	# Physically this means that this Bayes factor value for this points number
	# is nowhere achieved (and all available values are smaller)
	if zeta_t_lims < 0:
		print("Warning: zeta_t limits do not exist for these points number and K_L")
		return [[float('nan'), float('nan')], False]

	zeta_t_lims = np.sqrt(n_0 ** 2.0 / n / n_pi * zeta_t_lims) * np.asarray([-1, 1])
	print (zeta_t_lims)
	return [zeta_t_lims, True]
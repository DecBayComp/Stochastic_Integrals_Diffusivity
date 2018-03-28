

import numpy as np
from scipy import integrate


def calculate_marginalized_integral(zeta_t, zeta_sp, n, n_pi):
	"""
	Calculate the lambda integral for the given values of zeta_t and zeta_sp. 
	They must be 1D lists of same length.
	n and n0 are scalars.
	"""

	# Constants
	n0 = n + n_pi
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

		# If input is 0, no need to actually integrate
		if np.isclose(zeta_sp_cur, 0, rtol = tol, atol = tol):
			result = (1.0 + n / n0 * zeta_t_cur ** 2.0) ** ((2.0 - n0) / 2.0)
			results.append([result, 0])
			continue

		# Check if a break-point needs to be added
		break_lambda = zeta_t_cur/zeta_sp_cur
		if break_lambda < 0.0 or break_lambda > 1.0:
			break_lambda = None
		else:
			break_lambda = [break_lambda]

		# Define the integrand function
		def integrate_me(l):
			return (1.0 + n / n0 * (l * zeta_sp_cur - zeta_t_cur) ** 2.0) ** ((2.0 - n0) / 2.0)

		# Perform integration
		result = integrate.quad(integrate_me, 0.0, 1.0, points = break_lambda, full_output = 0)

		results.append(result)

	return(results)







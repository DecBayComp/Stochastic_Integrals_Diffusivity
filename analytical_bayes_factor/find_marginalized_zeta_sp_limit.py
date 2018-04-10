

from calculate_marginalized_integral import calculate_marginalized_integral
from log_C import log_C as log_C_func
import numpy as np
import scipy.optimize as optimize


def find_marginalized_zeta_sp_limit(n, n_pi, K):
	"""
	Find the max. abs. value of zeta_sp that provides the given value of the Bayes factor.
	It must lie on the line zeta_t = zeta_sp / 2.
	"""
	n0 = n + n_pi
	log_C = log_C_func(n0, n_pi)
	factor = 2.0
	max_attempts = 50

	# Construct the function to be optimized
	def solve_me(zeta_sp):
		# res = calculate_marginalized_integral([zeta_sp / 2.0], [zeta_sp], n, n_pi)
		res = (K - np.exp(log_C) / calculate_marginalized_integral([zeta_sp / 2.0], [zeta_sp], n, n_pi)[0][0])
		return res

	# Find the initial search interval
	# The left boundary zeta_sp = 0 is clear giving the min value of the integral the min value of the solve_me function
	# It's hard to define scale for the right border, but it shouldn't be a problem, so let's start with one
	search_interval = np.asarray([0, 1])
	min_value = K - np.exp(log_C)

	# Check that the min. value is positive
	if min_value < 0:
		# print ("Warning: zeta_sp limits do not exist for the given Bayes factor K = %.3g and n = %i" % (K, n))
		return 0

	# Start searching for the initial search interval
	attempt = 1
	bl_interval_found = False
	while attempt <= max_attempts:
		if solve_me(factor ** attempt) * min_value < 0:
			bl_interval_found = True
			search_interval[1] = factor ** attempt
			break
		attempt += 1

	if not bl_interval_found:
		raise RuntimeError("Unable to find initial search interval within %i attempts" % (max_attempts))

	# Find the actual root
	res = optimize.brentq(solve_me, search_interval[0], search_interval[1])
	return res

















# Copyright © 2018, Alexander Serov


from calculate_marginalized_integral import calculate_marginalized_integral
import numpy as np
import scipy


def calculate_minimal_n(zeta_t, zeta_sp, V, V_pi):
	"""
	Calculate the minimal number of jumps per bin to obtain strong evidence for a conservative 
	or a spurious force in the Bayes factor.

	Input:
	zeta_t, zeta_sp --- vectors of length 2 (x, y) in one bin

	Output:
	min_n --- minimal number of jumps to obtain strong evidence for one or another model.
	"""

	# Local constants
	n_pi = 4
	dim = 2
	B_threshold = 10.0
	increase_factor = 2
	max_attempts = 20
	xtol = 0.1
	rtol = 0.1


	# Initialize
	lg_B_threshold = np.log10(B_threshold)
	u = V_pi / V
	def eta(n):
		return np.sqrt(n_pi / (n + n_pi))

	def p(n):
		return dim * (n + n_pi + 1) / 2 - 2

	def v0(n):
		return 1.0 + n_pi / n * u

	# Define the Bayes factor
	def lg_B(n):
		upstairs = calculate_marginalized_integral(zeta_t = zeta_t, zeta_sp = zeta_sp, p = p(n),
				v = v0(n), E = eta(n)**2.0)
		downstairs = calculate_marginalized_integral(zeta_t = zeta_t, zeta_sp = zeta_sp, p = p(n),
			v = v0(n), E = 1.0)
		lg_B = dim * np.log10(eta(n)) + np.log10(upstairs) - np.log10(downstairs)
		return lg_B


	# Find initial search interval
	bl_found = False
	for attempt in range(max_attempts):
		n = increase_factor ** attempt
		if np.abs(lg_B(n)) >= lg_B_threshold:
			bl_found = True
			break

	if not bl_found:
		return np.nan

	if attempt == 0:
		return 1

	# Construct the initial interval
	n_interval = increase_factor ** np.asarray([0.0, attempt])
	lg_B_interval = [lg_B(n) for n in n_interval]
	
	# Identify whether the spurious or the conservative force model is favored by this n
	sign = np.sign(lg_B_interval[1])

	# Define the function being optimized
	def solve_me(n):
		return lg_B(n) - sign * lg_B_threshold

	# Find root
	min_n = scipy.optimize.brentq(solve_me, n_interval[0], n_interval[1])

	# Round
	min_n = np.ceil(min_n)
	
	return min_n

	








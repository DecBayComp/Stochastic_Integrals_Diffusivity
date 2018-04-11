

from calculate_marginalized_integral import calculate_marginalized_integral
import numpy as np

def calculate_bayes_factors(zeta_ts, zeta_sps, ns, Vs, Vs_pi):
	"""
	Calculate the marginalized Bayes factor for the presence of the conservative force in a set of bins.
	
	Input:
	zeta_ts --- signal-to-noise ratios for the total force = dx_mean / sqrt(var(dx)) in bins. Size: M x 2,
	zeta_sps --- signal-to-noise ratios for the spurious force = grad(D) / sqrt(var(dx)) in bins. Size: M x 2,
	ns --- number of jumps in each bin. Size: M x 1,
	Vs --- jump variance in bins = E((dx - dx_mean) ** 2). Size: M x 1,
	Vs_pi --- jump variance in all other bins relative to the current bin. Size: M x 1.

	Output:
	Bs, bl_forces

	Bs --- Bayes factor values in the bins. Size: M x 1,
	bl_forces --- Booleans for the presence of a conservative force. Bayes factors are thresholded at a level corresponding to strong evidence. Size: M x 1.

	Notation:
	M --- number of bins.

	"""

	# Constants
	n_pi = 4
	dim = 2
	B_threshold = 10	# corresponds to strong evidence for the conservative force

	# Convert input to numpy
	zeta_ts = np.asarray(zeta_ts)
	zeta_sps = np.asarray(zeta_sps)
	Vs = np.asarray(Vs)
	Vs_pi = np.asarray(Vs_pi)
	ns = np.asarray(ns)

	# Check that the 2nd dimension has size 2
	if np.shape(zeta_ts)[1] != 2 or np.shape(zeta_sps)[1] != 2:
		raise VsalueError("zeta_ts and zeta_sps must be matrices of size (M x 2)")

	# Initialize
	M = len(ns)
	etas = np.sqrt(n_pi / (ns + n_pi))
	us = np.divide(Vs_pi, Vs)
	v0s = 1.0 + np.multiply(n_pi / ns, us)
	pows = dim * (ns + n_pi + 1.0) / 2.0 - 2.0

	# Calculate 
	lg_Bs = np.zeros((M, 1))
	for i in range(M):
		upstairs = calculate_marginalized_integral(zeta_t = zeta_ts[i, :], zeta_sp = zeta_sps[i, :], p = pows[i],
			v = v0s[i], E = etas[i]**2.0)
		downstairs = calculate_marginalized_integral(zeta_t = zeta_ts[i, :], zeta_sp = zeta_sps[i, :], p = pows[i],
			v = v0s[i], E = 1.0)
		lg_Bs[i] = dim * np.log10(etas[i]) + np.log10(upstairs) - np.log10(downstairs)
	
	# Threshold
	bl_force = lg_Bs >= np.log10(B_threshold)

	return (10.0 ** lg_Bs, bl_force)





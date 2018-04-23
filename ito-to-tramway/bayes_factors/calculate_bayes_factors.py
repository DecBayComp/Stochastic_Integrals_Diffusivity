# Copyright © 2018, Alexander Serov


from .calculate_marginalized_integral import calculate_marginalized_integral
from .calculate_minimal_n import calculate_minimal_n
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
	Bs, forcess

	Bs --- Bayes factor values in the bins. Size: M x 1,
	forces --- Returns 1 if there is strong evidence for the presence of a conservative forces, 
	-1 for strong evidence for 	a spurious force, and 0 if the is not enough evidence. Size: M x 1.

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
	lg_Bs = np.zeros((M, 1)) * np.nan
	min_ns = np.zeros((M, 1), dtype = int) * np.nan
	for i in range(M):
		try:
			upstairs = calculate_marginalized_integral(zeta_t = zeta_ts[i, :], zeta_sp = zeta_sps[i, :], p = pows[i],
				v = v0s[i], E = etas[i]**2.0)
			downstairs = calculate_marginalized_integral(zeta_t = zeta_ts[i, :], zeta_sp = zeta_sps[i, :], p = pows[i],
				v = v0s[i], E = 1.0)
			lg_Bs[i] = dim * np.log10(etas[i]) + np.log10(upstairs) - np.log10(downstairs)
			min_ns[i] = calculate_minimal_n(zeta_ts[i, :], zeta_sp = zeta_sps[i, :], n0 = ns[i], V = Vs[i], 
				V_pi = Vs_pi[i])
		except:
			print("Warning: Detected data error in bin %i. Skipping bin." % i)
	
	# Threshold into 3 categories: strong evidence for either of the models and insufficient evidence
	forces = 1 * (lg_Bs	>= np.log10(B_threshold)) - 1 * (lg_Bs <= -np.log10(B_threshold))
	# forces = np.zeros_like(lg_Bs, dtype = int)
	# forces [lg_Bs	>= np.log10(B_threshold)] = 1
	# forces [lg_Bs	<= - np.log10(B_threshold)] = -1

	return (10.0 ** lg_Bs, forces, min_ns)





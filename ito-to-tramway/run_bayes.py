

from calculate_bayes_factors import calculate_bayes_factors
from calculate_minimal_n import calculate_minimal_n
import numpy as np

# Initialize data
zeta_sps = [[1, 1], [2, 3]]
zeta_ts = [[2, 1.54], [-4, 0]]
ns = [10, 3]
Vs = [1, 0.2]
Vs_pi = [0.3, 0.5]


# Run calculations
# calculate_bayes_factors(zeta_ts = zeta_ts, zeta_sps = zeta_sps, ns = ns, Vs = Vs, Vs_pi = Vs_pi)
calculate_minimal_n(zeta_t = np.asarray(zeta_ts)[0, :], zeta_sp = np.asarray(zeta_sps)[0, :], 
	V = Vs[0], V_pi = Vs_pi [0])

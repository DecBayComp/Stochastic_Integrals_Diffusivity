

from calculate_zeta_t_lims import calculate_zeta_t_lims
import numpy as np
import sys



## Constants
lambdas = (0.5, 1.0)
zeta_t_steps = 20
n_pi = 4
# K_L_level = 1 # not used
n = 1000


# Calculate zeta_t border for lambda > 0
n_0 = n + n_pi
[zeta_t_lims, bl_lims_exist] = calculate_zeta_t_lims(n_0, n, n_pi)
if not bl_lims_exist:
	print("Calculations aborted!")
	sys.exit(0)

# Linear space over zeta_t
zeta_t_points = np.linspace(zeta_t_lims[0], zeta_t_lims[1], zeta_t_steps)

# # Initialize root arrays
lambdas_count = len(lambdas)
klf_roots_size = (zeta_t_steps, lambdas_count, 2)
# klf_roots = np.zeros(() , dtype = float)










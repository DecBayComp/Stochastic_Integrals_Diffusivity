

from calculate_zeta_t_lims import calculate_zeta_t_lims
from calculate_zeta_sp_roots import calculate_zeta_sp_roots
import numpy as np
from plot_K_L import plot_K_L
import sys



## Constants
lambs = (0.5, 1.0)
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

# Initialize root arrays
lambs_count = len(lambs)
zeta_sp_roots_size = (zeta_t_steps, lambs_count, 2)
zeta_sp_roots = np.zeros(zeta_sp_roots_size, dtype = np.float)

# Calculate roots
for lamb_ind in range(lambs_count):
	lamb = lambs[lamb_ind]
	
	# Not on zeta_t borders
	for zeta_ind in range(zeta_t_steps):

		zeta_t = zeta_t_points[zeta_ind]
		zeta_sp_roots[zeta_ind, lamb_ind, :] = calculate_zeta_sp_roots(zeta_t, lamb, n_0, n, n_pi)
	
	# On the borders the formula is different
	zeta_sp_roots[0, lamb_ind, :] = zeta_t_points[0] / lamb
	zeta_sp_roots[-1, lamb_ind, :] = zeta_t_points[-1] / lamb

print(zeta_sp_roots)











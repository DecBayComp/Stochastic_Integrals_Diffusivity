

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def plot_zeta_sp_borders(ns, zeta_sp_borders):

	ns_count = len(ns)
	Ks_count = np.size(zeta_sp_borders, 0)
	plt.clf()

	for K_ind in range(Ks_count):
		plt.plot(ns, zeta_sp_borders[K_ind])

	plt.xlabel("n")
	plt.ylabel("max(zeta_sp)")
	plt.show()






from constants import *
import matplotlib.pyplot as plt
import numpy as np


def plot_K_L(zeta_t_points_abs, zeta_sp_roots):
	"""
	Make the parametric plot of the local Bayes factor
	"""
	lambs_count = np.size(zeta_sp_roots, 0)
	# branch_ind = 0
	# fig = mpl.pyplot.figure()

	# Make a cyclic plot mesh
	zeta_t_points_back = np.flipud(zeta_t_points_abs[0:-1])
	print(zeta_t_points_abs)
	print(zeta_t_points_back)
	zeta_t_plot_mesh = np.append(zeta_t_points_abs, zeta_t_points_back)

	plt.clf()

	lamb_ind = 0
	for lamb_ind in range(lambs_count):
	# 	for branch_ind in [0, 1]:
		# y_axis = np.asarray(zeta_t_plot_mesh) / np.asarray(zeta_sp_roots[lamb_ind, :])
		# y_axis = 1 / y_axis
		y_axis =  np.asarray(zeta_sp_roots[lamb_ind, :])

		# Positive zeta_t branch
		plt.plot(zeta_t_plot_mesh, y_axis, color = color_sequence[lamb_ind + 1])

		# Negative zeta_t branch
		plt.plot(- zeta_t_plot_mesh, y_axis, color = color_sequence[lamb_ind + 1])

		# branch_ind = 1
		# plt.plot(zeta_t_points, zeta_sp_roots[:, lamb_ind, branch_ind], color = color_sequence[lamb_ind + 1])
	
	# Adjust
	# plt.axis([zeta_t_points[0] * 1.05, zeta_t_points[-1] * 1.05, -2, 2])
	plt.xlabel('zeta_t')
	plt.ylabel('zeta_sp')

	plt.show()
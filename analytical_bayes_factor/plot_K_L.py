

from constants import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def plot_K_L(zeta_sps, zeta_t_roots, ns):
	"""
	Make the parametric plot of the local Bayes factor
	"""
	# Constants
	alpha = 0.18
	font_size = 8
	pagewidth_in = 6.85
	dpi = 120
	figsize = np.asarray([1, 0.4]) * pagewidth_in # in inches

	lambs_count = np.size(zeta_t_roots, 1)
	ns_count = len(ns)
	# branch_ind = 0
	# fig = mpl.pyplot.figure()

	# # Make a cyclic plot mesh
	# zeta_t_points_back = np.flipud(zeta_t_points_abs[0:-1])
	# print(zeta_t_points_abs)
	# print(zeta_t_points_back)
	# zeta_t_plot_mesh = np.append(zeta_t_points_abs, zeta_t_points_back)

	# fig = plt.figure(1)
	# plt.clf()
	matplotlib.rcParams.update({'font.size': font_size})
	plt.clf()
	plt.rc('text', usetex = True)
	# plt.style.use("axes.labelsize = font_size")

	fig, axarr = plt.subplots(1, 3, num = 1, sharey = True)

	# Adjust
	fig.set_dpi(dpi)
	fig.set_figwidth(figsize[0])
	fig.set_figheight(figsize[1])
	

	# plt.subplots_adjust(wspace=0, hspace=0)
	# ax.labelsize = font_size


	for n_ind in range(ns_count):
		ax = axarr[n_ind]
		n = ns[n_ind]

		# Plot fixed-lambda
		for lamb_ind in range(lambs_count - 1):
		# 	for branch_ind in [0, 1]:
			# y_axis = np.asarray(zeta_t_plot_mesh) / np.asarray(zeta_sps[lamb_ind, :])
			# y_axis = 1 / y_axis
			# y_axis =  np.asarray(zeta_sps[lamb_ind, :])

			# Positive zeta_t branch
			# branch_ind = 0
			# y_axis = zeta_t_roots[lamb_ind, branch_ind, :]
			print(zeta_t_roots[n_ind, lamb_ind, 0, :])
			ax.fill_between(zeta_sps, zeta_t_roots[n_ind, lamb_ind, 0, :], zeta_t_roots[n_ind, lamb_ind, 1, :],
				color = color_sequence[lamb_ind], alpha = alpha)

			# Negative zeta_t branch
			# branch_ind = 1
			# y_axis = zeta_t_roots[lamb_ind, branch_ind, :]
			# plt.plot(zeta_sps, y_axis, color = color_sequence[lamb_ind])

			# branch_ind = 1
			# plt.plot(zeta_t_points, zeta_sps[:, lamb_ind, branch_ind], color = color_sequence[lamb_ind + 1])


		lamb_ind = lambs_count - 1
		# 	for branch_ind in [0, 1]:
		# y_axis = np.asarray(zeta_t_plot_mesh) / np.asarray(zeta_sps[lamb_ind, :])
		# y_axis = 1 / y_axis
		# y_axis =  np.asarray(zeta_sps[lamb_ind, :])

		# Positive zeta_t branch
		branch_ind = 0
		y_axis = zeta_t_roots[n_ind, lamb_ind, branch_ind, :]
		ax.plot(zeta_sps, y_axis, color = color_sequence[lamb_ind])

		# Negative zeta_t branch
		branch_ind = 1
		y_axis = zeta_t_roots[n_ind, lamb_ind, branch_ind, :]
		ax.plot(zeta_sps, y_axis, color = color_sequence[lamb_ind])

		# branch_ind = 1
		# plt.plot(zeta_t_points, zeta_sps[:, lamb_ind, branch_ind], color = color_sequence[lamb_ind + 1])
		
		# Adjust
		# plt.axis([zeta_t_points[0] * 1.05, zeta_t_points[-1] * 1.05, -2, 2])
		ax.set_xlabel('$\zeta_\mathrm{sp}$')
		
		ax.set_title("n = %i" % (n))


	axarr[0].set_ylabel('$\zeta_\mathrm{t}$')
	plt.tight_layout()
	# plt.subplots_adjust(wspace=0, hspace=0, left = 0.0, right = 1.0)
	plt.show()
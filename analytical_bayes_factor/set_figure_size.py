

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def set_figure_size(num, rows, page_width_frac):
	pagewidth_in = 6.85
	font_size = 8
	dpi = 120

	figsize = np.asarray([page_width_frac, 0.33 * rows]) * pagewidth_in # in inches
	
	# Set default font size and LaTeX usage
	matplotlib.rcParams.update({'font.size': font_size})
	plt.rc('text', usetex = True)

	# Create and return figure handle
	fig = plt.figure(num)
	fig.clf()

	# Set figure size and dpi
	fig.set_dpi(dpi)
	fig.set_figwidth(figsize[0])
	fig.set_figheight(figsize[1])

	# fig.tight_layout()

	# Return figure handle
	return (fig)
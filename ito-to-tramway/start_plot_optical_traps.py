"""
Calculate and plot the Bayes factors in optical traps provided by Max. The idea is to keep the number of jumps per bin the same, but change the power of the optical trap. Due to this, there will be different level of evidence for the active forces in the different cases, although in all cases there is confinement. In the weakest cases, the confienement might be due to another stochastic integral interpretation.
"""

try:
    has_run
except NameError:
    %matplotlib
    %load_ext autoreload
    %autoreload 2

    has_run = 1
else:
    print("Graphic interface NOT re-initialized")

from calculate import calculate
from constants import optical_traps_data_folder, optical_data_sets, optical_traps_points_per_bin, MACHINE_PRECISION, optical_traps_dt, localization_error, optical_power_mW
from bayes_comparison_optical_traps.load_trajectory import load_trajectory
import matplotlib
import matplotlib.pyplot as plt
import os
import pandas as pd
from set_figure_size import set_figure_size
from snr import infer_snr
from tqdm import trange
from tramway.helper import *
from tramway.plot.map import scalar_map_2d

# def plot_optical_traps():

# %% Constants
font_size = 8
alpha = 1.0
linewidth = 0.1
dpi = 100


# %%
cells_combined = []
log10_Bs_combined = []

for file_ind in trange(len(optical_data_sets)):
    file = optical_data_sets[file_ind]

    # Load current trajectory
    filename = 'Brownien_%s' % (file)
    trajectory, rwa_fullpath = load_trajectory(filename)

    # Tessellate with a fixed number of jumps per bin
    mesh = tessellate(trajectory, 'kmeans', label='equal', verbose=False,
                      avg_location_count=optical_traps_points_per_bin, ref_distance=MACHINE_PRECISION, prune=0, force=True)

    # Perform some modifications
    rwa_data = Analyses(trajectory)
    rwa_data.add(mesh, label='equal')
    equal_tree = rwa_data['equal']
    equal_tree_distr = distributed(equal_tree.data)

    # Print the number of jumps per bin
    ns = np.asarray([len(equal_tree_distr[i]) for i in equal_tree_distr])
    # print(ns)

    # fig = plt.figure(1)
    # plt.clf()
    # cell_plot(mesh, num=1, figsize=None)

    # Perform inference and calculate the Bayes factors in each bin
    # -> May need to check the localization error value
    inference_result = equal_tree_distr.run(
        infer_snr, max_iter=50, localization_error=localization_error)
    equal_tree.add(Maps(inference_result, mode='snr'), label='snr')
    save_rwa(rwa_fullpath, rwa_data, force=True)

    # Calculate the Bayes factors
    output_folder = r".\temp"
    calculate(rwa_fullpath, output_folder, bl_produce_maps=False, dt=optical_traps_dt,
              snr_label='snr', localization_error=localization_error)

    # Get the calculated Bayes factors
    rwa_data = load_rwa(rwa_fullpath)
    cells = rwa_data['equal'].data
    n = rwa_data['equal']['snr'].data['n']
    log10_Bs = rwa_data['equal']['snr']['bayes_factors'].data.log10_B
    # Rename
    # log10_Bs.columns = [file_ind]

    # Save to structures
    log10_Bs_combined.append(log10_Bs)
    cells_combined.append(cells)


# %% Initialize the figure
# Initialize a figure to incorporate all bayes factor plots
num = 1
rows = 1
cols = 3
page_width_frac = 1.0
pagewidth_in = 6.85
height_factor = 0.9

fig, figsize = set_figure_size(
    num=num, rows=rows, page_width_frac=page_width_frac, height_factor=height_factor)
fig.clf()
_, axes = plt.subplots(nrows=rows, ncols=cols, num=num, sharex=True, sharey=True)

# Determine min and max values
vmin = np.nanmin([np.min(vals) for vals in log10_Bs_combined])
vmax = np.nanmax([np.max(vals) for vals in log10_Bs_combined])


# % Plot a map of the Bayes factors
for file_ind in range(len(optical_data_sets)):
    ax = axes[file_ind]
    plt.sca(ax)
    cells = cells_combined[file_ind]
    log10_Bs = log10_Bs_combined[file_ind]
    bl_colorbar = None
    if file_ind == len(optical_data_sets) - 1:
        bl_colorbar = 'nice'

    scalar_map_2d(cells=cells, values=log10_Bs, axes=ax, colormap='inferno',
                  alpha=alpha, linewidth=linewidth, aspect=None, clim=[vmin, vmax], colorbar=bl_colorbar)

    str_title = '%i mW' % (optical_power_mW[file_ind])
    ax.set_title(str_title)
    ax.set_xlabel('$x$, $\mu$m')
    if file_ind == 0:
        ax.set_ylabel('$y$, $\mu$m')

    # add label
    label_location = [0.025, 1.03]
    str_label = chr(ord('a') + file_ind)
    # plot_me.count += 1
    ax = fig.gca()
    ax.text(label_location[0], label_location[1],
            str_label, transform=ax.transAxes, fontsize=font_size)

    # Colorbar legend
    colorbar_legend = '$\log_{10} (B)$'
    if file_ind == len(optical_data_sets) - 1:
        ax = fig.axes[len(optical_data_sets)].set_ylabel(colorbar_legend, rotation=90)

fig.tight_layout()
fig.show()

# Save
fig_basename = "optical_traps_bayes_factors"
fig.savefig(fig_basename + ".pdf")
fig.savefig(fig_basename + ".png")

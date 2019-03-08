
"""
Calculate and plot the Bayes factors in optical traps provided by Max. The idea is to keep the number of jumps per bin the same, but change the intensity of the optical trap. Due to this, there will be different level of evidence for the active forces in the different cases, although in all cases there is confinement. In the weakest cases, the confienement might be due to another stochastic integral interpretation.
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

import os
import time

# from rwa import HDF5Store
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import trange

from bayes_comparison_optical_traps.load_trajectory import load_trajectory
from calculate import calculate
from constants import (MACHINE_PRECISION, optical_data_sets, optical_power_mW,
                       optical_traps_data_folder, optical_traps_dt,
                       optical_traps_points_per_bin, pagewidth_in)
from get_exterior_cells import get_exterior_cells
from set_figure_size import set_figure_size
from snr import infer_snr
from stopwatch import stopwatch
from tramway.helper import *
from tramway.plot.map import scalar_map_2d
from tramway.plot.mesh import *

# %% Constants
sigma = 1.0e-8  # um
font_size = 8
alpha = 1.0
linewidth = 0.1
dpi = 100
tessellation_alg = 'kohonen'
bl_recalculate_mesh = False
fig_basename = "optical_traps_bayes_factors"

label = "%s-%i" % (tessellation_alg, optical_traps_points_per_bin)

print("Using folder '%s'" % (optical_traps_data_folder))
print("Using label '%s'" % (label))


# %% Tessellate and calculate the Bayes factors
cells_combined = []
log10_Bs_combined = []
ns_combined = []

for file_ind in trange(len(optical_data_sets)):
    # file_ind = 0
    file = optical_data_sets[file_ind]

    # Load current trajectory
    filename = 'Brownien_%s' % (file)
    trajectory, rwa_fullpath = load_trajectory(filename)

    # Load .rwa and check for a pre-calculated mesh
    rwa_data = []
    if os.path.exists(rwa_fullpath):
        rwa_data = load_rwa(rwa_fullpath, verbose=False)

    if label not in rwa_data or bl_recalculate_mesh:
        # Tessellate with an approximately fixed number of jumps per bin
        t0 = time.time()
        mesh = tessellate(trajectory, tessellation_alg, label=label, verbose=True,
                          avg_location_count=optical_traps_points_per_bin, rel_max_distance=.5, knn=(optical_traps_points_per_bin - 10, optical_traps_points_per_bin + 10))
        print("Tessellation completed in % .1f s" % (time.time() - t0))

        # Save the new mesh into the .rwa file
        rwa_data = Analyses(trajectory)
        rwa_data.add(mesh, label=label)
        save_rwa(rwa_fullpath, rwa_data, force=True)
    else:
        print("Loaded a previously calculated tessellation mesh")

    equal_tree = rwa_data[label]
    equal_tree_distr = distributed(equal_tree.data)
    cells = rwa_data[label].data
    cells_len = cells.location_count.size

    # Calculate the Bayes factors
    with stopwatch('Bayes factor calculations'):
        calculate(rwa_fullpath, results_folder=None, bl_produce_maps=False,
                  snr_label='snr', sigma=sigma)

    # Get the calculated Bayes factors
    rwa_data = load_rwa(rwa_fullpath)
    log10_Bs = rwa_data[label]['snr']['bayes_factor'].data['lg_B']['lg_B'].copy()

    # Mark surface cells, where the results may not be accurate
    surface_cells = get_exterior_cells(cells)
    cells.tessellation.cell_label = np.ones(cells_len, dtype=bool)
    cells.tessellation.cell_label[surface_cells] = False
    log10_Bs[surface_cells] = np.nan

    # Save to structures
    log10_Bs_combined.append(log10_Bs)
    cells_combined.append(cells)


# %% Plot the Bayes factors together
num = 1
rows = 1
cols = 1
page_width_frac = 0.5
height_factor = 0.75
aspect = 'equal'
xlim = np.asarray([-1, 1]) * 9
ylim = xlim

fig, figsize = set_figure_size(
    num=num, rows=rows, page_width_frac=page_width_frac, height_factor=height_factor)

_, ax = plt.subplots(nrows=rows, ncols=cols, num=num, aspect=aspect)

# Determine the colorbar scale
min_lg_B = np.nanmin([np.nanmin(vals) for vals in log10_Bs_combined])
max_lg_B = 3  # np.nanmax([np.max(vals) for vals in log10_Bs_combined])

# % Plot a map of the Bayes factors
row = 0
for file_ind in trange(len(optical_data_sets)):
    ax.clear()

    cells = cells_combined[file_ind]
    log10_Bs = log10_Bs_combined[file_ind]
    bl_colorbar = 'nice'

    scalar_map_2d(cells=cells, values=log10_Bs, axes=ax, colormap='inferno',
                  alpha=alpha, linewidth=linewidth, aspect=aspect, clim=[min_lg_B, max_lg_B], colorbar=bl_colorbar)

    # Adjust
    str_title = '%i mW' % (optical_power_mW[file_ind])
    ax.set_title(str_title)

    ax.set_xlabel('$x$, $\mu$m')
    ax.set_ylabel('$y$, $\mu$m')

    plt.xlim(xlim)
    plt.ylim(ylim)
    ticks = [-5, 0, 5]
    plt.xticks(ticks)
    plt.yticks(ticks)

    # # Add letter label
    # label_location = [0.025, 1.03]
    # str_label = chr(ord('a') + file_ind)
    # ax.text(label_location[0], label_location[1], str_label,
    #         transform=ax.transAxes, fontsize=font_size)

    # Colorbar title
    colorbar_title = '$\log_{10} K$'
    fig.axes[rows * cols].set_ylabel(colorbar_title, rotation=90)

    fig.tight_layout()
    fig.show()

    # Save
    figname = f'{fig_basename}_{optical_power_mW[file_ind]:d}'
    fig.savefig(figname + ".pdf", bbox_inches='tight', pad_inches=0)
    fig.savefig(figname + ".png", bbox_inches='tight', pad_inches=0)

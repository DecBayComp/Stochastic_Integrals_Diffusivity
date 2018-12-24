try:
    bl_has_run

except Exception:
    %load_ext autoreload
    %autoreload 2
    bl_has_run = True

import datetime
import os.path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from tramway.helper import (Analyses, infer, load_rwa, load_xyt, map_plot,
                            save_rwa, tessellate)
from tramway.inference import bayes_factors

# %matplotlib inline

# files and parameters
input_folder = './input'

# trajectory_file = 'Brownien_P1.txt'
# trajectory_file = 'Brownien_P4.txt'
# trajectory_file = 'VLP_WT_2_2.txt'
trajectory_file = 'VLP_WT_2_full.txt'

localization_error = 0
localization_error_BF = localization_error**2

cell_radii = []  # 0.8, 0.4]  # .03]  # , .09, .27, .81]
location_counts = [20]  # 1000, 500]  # 20, 80, 320, 1280]
B_thresholds = [10, 2]

trajectory_file = os.path.join(input_folder, trajectory_file)
rwa_file = os.path.splitext(trajectory_file)[0] + '.rwa'


# %%
if not os.path.isfile(rwa_file):
    # load the trajectories
    if 'Brownien' in trajectory_file:
        # time step (please check)
        dt = 1.0 / 65536
        # read the table
        xyt = load_xyt(trajectory_file, columns=['x', 'y', 'H'])
        # add missing columns
        xyt['n'] = np.ones(xyt.shape[0])
        xyt['t'] = np.arange(dt, (xyt.shape[0] + 1) * dt, dt)
    elif 'VLP' in trajectory_file:
        xyt = load_xyt(trajectory_file, columns=['n', 'x', 'y', 't'])
    else:
        xyt = load_xyt(trajectory_file)
    analysis_tree = Analyses(xyt)
    save_rwa(rwa_file, analysis_tree, force=True)
else:
    analysis_tree = load_rwa(rwa_file)

# tessellate
for radius in tqdm(cell_radii):
    label = 'hexagon_r_{:.2f}'.format(radius)
    if label not in analysis_tree:
        cells = tessellate(analysis_tree, 'hexagon', ref_distance=radius,
                           min_location_count=0, label=label)
        save_rwa(rwa_file, analysis_tree, force=True)

for count in tqdm(location_counts):
    label = 'kmeans_n_{:d}'.format(count)
    if label not in analysis_tree:
        cells = tessellate(analysis_tree, 'kmeans', ref_distance=0, avg_location_count=count,
                           knn=(round(.75 * count), round(1.25 * count)), prune=False, label=label)
        save_rwa(rwa_file, analysis_tree, force=True)

# save to file
# save_rwa(rwa_file, analysis_tree, force=True)

# %%
snr_label = 'snr'
bayes_factor_label = 'bayes_factor'

analysis_tree = load_rwa(rwa_file)

for mesh in analysis_tree.labels:
    for B_threshold in B_thresholds:
        print("Bayes factors for mesh '{mesh}' with B threshold = {B_threshold}".format(
            mesh=mesh, B_threshold=B_threshold))

        # infer forces and diffusivity
        if snr_label not in analysis_tree[mesh]:
            infer(analysis_tree, 'snr', input_label=mesh, output_label=snr_label,
                  max_iter=50, localization_error=localization_error)
            save_rwa(rwa_file, analysis_tree, force=True)

        output_label = bayes_factor_label + str(B_threshold)
        if output_label not in analysis_tree[mesh][snr_label]:
            BF_results = infer(analysis_tree, 'bayes_factor',                            input_label=[mesh, snr_label],
                               output_label=output_label, localization_error=localization_error_BF, B_threshold=B_threshold)
            save_rwa(rwa_file, analysis_tree, force=True)
        else:
            BF_results = analysis_tree[mesh][snr_label][output_label].data

        lg_Bs = BF_results['lg_B']
        forces = BF_results['force']

        cells = analysis_tree[mesh].data
        map_plot(lg_Bs, cells=cells,
                 show=True, clip=False, colormap='inferno', alpha=1, linewidth=0.1, figsize=(6.85, 2.29), dpi=100, aspect='equal', colorbar='nice')

        map_plot(forces, cells=cells,
                 show=True, clip=False, colormap='inferno', alpha=1, linewidth=0.1, figsize=(6.85, 2.29), dpi=100, aspect='equal', colorbar='nice')

    # fig = plt.gcf()
    # fig.savefig(os.path.join('input', 'out-lg_B.png'))

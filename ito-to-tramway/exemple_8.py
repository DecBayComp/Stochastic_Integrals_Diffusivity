
%matplotlib inline

# %%
# files and parameters

trajectory_file = 'Brownien_P1.txt'

cell_radii = [.03]  # , .09, .27, .81]
location_counts = []  # 20]  # , 80, 320, 1280]

#

import os.path

import numpy as np

from tramway.helper import *
from tramway.inference.bayes_factors import calculate_bayes_factors

trajectory_file = os.path.join('input', trajectory_file)
rwa_file = os.path.splitext(trajectory_file)[0] + '.rwa'


# %%
# load the trajectories

if 'Brownien' in trajectory_file:

    # time step (please check)
    dt = 1.0 / 65536
    # read the table
    xyt = load_xyt(trajectory_file, columns=['x', 'y', 'H'])
    # add missing columns
    xyt['n'] = np.ones(xyt.shape[0])
    xyt['t'] = np.arange(dt, (xyt.shape[0] + 1) * dt, dt)

else:

    xyt = load_xyt(trajectory_file)

# %%
# tessellate

analysis_tree = Analyses(xyt)

mesh_labels = []

for radius in cell_radii:
    cells = tessellate(xyt, 'hexagon', ref_distance=radius, min_location_count=0)
    analysis_tree.add(cells, label='hexagon_r_{:.2f}'.format(radius))

for count in location_counts:
    cells = tessellate(xyt, 'kmeans', ref_distance=0, avg_location_count=count,
                       knn=(round(.75 * count), round(1.25 * count)), prune=False)
    analysis_tree.add(cells, label='kmeans_n_{:d}'.format(count))

# %%
# save to file

save_rwa(rwa_file, analysis_tree, force=True)

# %%
# or start from here and load the .rwa file

analysis_tree = load_rwa(rwa_file)

# %%
# play

trajectories = analysis_tree.data

for mesh_label in analysis_tree.labels:

    cells = analysis_tree[mesh_label].data

    # infer forces and diffusivity
    if snr_label not in analysis_tree[mesh]:
        infer(analysis_tree, 'snr', input_label=mesh, output_label=snr_label,
              max_iter=50, localization_error=localization_error)
        save_rwa(rwa_file, analysis_tree, force=True)

    snr = analysis_tree[mesh][snr_label].data

    # calculate Bayes factor here
    Bs, forces, min_ns = calculate_bayes_factors(
        zeta_ts=zeta_ts, zeta_sps=zeta_sps, ns=ns, Vs=Vs, Vs_pi=Vs_prior)

    # save to the analysis tree

    # plot Bayes factor maps (and save them?)
    my_map = pd.DataFrame(np.log10(Bs), index=n.index, columns=['$\log_{10}(B)$'])
    map_plot(my_map, cells=cells,
             show=False, clip=False, colormap=colormap, alpha=alpha, linewidth=linewidth, figsize=figsize, dpi=dpi, aspect='equal', colorbar=colorbar)

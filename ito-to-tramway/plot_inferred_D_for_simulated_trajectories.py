
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

import numpy as np
import pandas as pd

from calculate import calculate
from constants import dt
from plot_D_profile import plot_D_profile
# import this
from tesselate_and_infer import tesselate_and_infer
from tramway.core.hdf5.store import load_rwa

# %% Produce a diffusivity map for one of the simulated trajectories for the article
filename = 'sim_data_000000353'
file = os.path.join(
    r'D:\Google Drive\git\Stochastic_Integrals_Diffusivity\ito-to-tramway\input\diffusivity_map_for_article_sim_trajectory', filename + '.csv')
output_folder = r'D:\Google Drive\git\Stochastic_Integrals_Diffusivity\ito-to-tramway\input\diffusivity_map_for_article_sim_trajectory\result'
dat_file = os.path.join(output_folder, filename + '_gwr.dat')

tesselate_and_infer(file, sigma2=1e-8, load=1)

# %%
# Output data and the diffusivity map figure
calculate(file, output_folder, bl_produce_maps=True,
          snr_label='snr', sigma=1e-8, ticks=True, page_width_frac=0.5)

# % Plot inferred 1D diffusivity profile along y0 line
y0 = 0.5
D_inferred_data = pd.read_csv(dat_file).loc[:, ['D', 'x_center', 'D_CI_lower', 'D_CI_upper']]
analysis_tree = load_rwa(file[:-4] + '.rwa')
get_coordinates = analysis_tree['gwr'].data.tessellation.cell_index

# Identify bins located on y=y0 line
x_mesh = np.linspace(0, 1, num=100)
probe_coordinates = [[x, y0] for x in x_mesh]
cells_on_line = set(get_coordinates(np.array(probe_coordinates)))
# Filter data
D_on_line = D_inferred_data.loc[cells_on_line, :]

plot_D_profile(D_on_line)

"""
This script collects .rwa files in subfolders on the network and creates an arguments file with their names.
The goal is to analyze J-B's simulated trajectories
"""

try:
    bl_has_run

except Exception:
    import matplotlib
    matplotlib.use('Agg')  # enable for console runs with no displays
    # %matplotlib
    %load_ext autoreload
    %autoreload 2
    bl_has_run = True

import copy
import glob
import logging
import os

# %%
import numpy as np
from scipy.special import gammaln
from tqdm import tqdm

from calculate import calculate
from load_trajectory import load_trajectory
from reinit_folder import reinit_folder
from tramway.core.analyses.lazy import Analyses
from tramway.helper import RWAStore

# % Constants
results_folder = "bayes_factors"
kwargs = {}


# root_path = r"\\157.99.40.171\@Dbc\LAB_shared_stuff\Francois_Laurent\tests_tramway\numerical_trajectories_no_box\snr"
# root_path = r"\\157.99.40.171\@Dbc\LAB_shared_stuff\Francois_Laurent\tests_tramway\numerical_trajectories"

# # Miscellaneous experiments
# root_path = r"\\157.99.40.171\@Dbc\LAB_shared_stuff\Francois_Laurent\tests_tramway\misc_experiments"
# dt = 0.05  # s
# snr_label = 'snr'
# localization_error = 0.03
# bl_recursive = False

# # optical tweezers
# root_path = r"\\atlas.pasteur.fr\@Dbc\LAB_shared_stuff\Francois_Laurent\tests_tramway\optical_tweezers"
# dt = 1 / 65636  # s
# snr_label = 'snr(mu=0)'
# localization_error = 0.0

# # optical tweezers local
# root_path = r'./input/optical_tweezers'
# snr_label = 'snr'
# localization_error = 0.0


# # Full viral capside
# # root_path = r"\\157.99.40.171\@Dbc\LAB_shared_stuff\Francois_Laurent\tests_tramway\misc_experiments\full"
# root_path = r'./input/vlp_full_francois'
# # root_path = r'./input/vlp_full'
# snr_label = 'snr'
# localization_error = 0.03**2


# VLP region
# root_path = r"\\atlas.pasteur.fr\@Dbc\LAB_shared_stuff\Francois_Laurent\tests_tramway\misc_experiments\vlp_2_2"
# root_path = r"./input/vlp_2_2"
# snr_label = 'snr'
# sigma = 0.03    # um


# Neuro-receptors
root_path = r".\input\neuro-receptors"
# root_path = r".\input\neuro-receptors-post"
# root_path = r".\input\neuro-receptors-pre-select"
# root_path = r".\input\neuro-receptors-post-select"
snr_label = 'snr'
sigma = 0.03    # um
# extension = '.txt'
extension = '.trxyt'
# xcenter, ycenter = 1.0, 1.0
# width = 0.5
kwargs = {
    'method': 'hexagon', 'avg_distance': 0.1, 'min_location_count': 0,
    # 'method': 'kohonen', 'min_probability': 20,  # 'avg_probability': 20,
    # 'avg_distance': 0.1,
    'reset_origin': False,
    # 'xlims': np.array([-1, 1]) / 2 * width + xcenter,
    # 'ylims': np.array([-1, 1]) / 2 * width + ycenter
}


# # optical tweezers local
# root_path = r'./input/optical_tweezers-test'
# snr_label = 'snr'
# sigma = 0.0
# dt = 1 / 65636  # s
# extension = '.txt'
# kwargs = {'method': 'hexagon', 'avg_distance': 2, 'min_location_count': 0}


# Varied beads lattice
root_path = r".\input\beads_lattice"
snr_label = 'snr'
sigma = 0.0    # um
extension = '.csv'
# xcenter, ycenter = 1.0, 1.0
# width = 0.5
kwargs = {
    'method': 'hexagon', 'avg_distance': 0.1, 'min_location_count': 0,
    # 'method': 'kohonen', 'min_probability': 20,  # 'avg_probability': 20,
    # 'avg_distance': 0.1,
    'reset_origin': False,
    # 'xlims': np.array([-1, 1]) / 2 * width + xcenter,
    # 'ylims': np.array([-1, 1]) / 2 * width + ycenter
}

# , 'rel_avg_distance': 0.01
# %%
bl_recursive = False
recalculate = True
if not extension:
    extension = '.txt'


# %% Parse all subdirectories to extract .rwa files and store

if bl_recursive:
    file_list = [f for f in glob.iglob(os.path.join(
        root_path, "/**/*" + extension), recursive=True)]
else:
    file_list = [f for f in glob.iglob(os.path.join(root_path, r"*" + extension), recursive=False)]

folder, _ = os.path.split(file_list[0])
results_path = os.path.join(folder, results_folder)
reinit_folder(results_path)

# %%
np.random.seed()
if len(file_list) == 0:
    logging.warning('Input file list empty')
else:
    for file in tqdm(file_list):
        folder, _ = os.path.split(file)
        print("Processing file: %s" % file)
        # analysis_tree = load_trajectory(file, txt_extension=extension)
        # print(analysis_tree)
        calculate(csv_file=file, results_folder=results_path, bl_produce_maps=True,
                  snr_label=snr_label, sigma=sigma, recalculate=recalculate, ticks=True, txt_extension=extension, pdf=True, **kwargs)

# #
# N = 100
# m = 50
# (gammaln(N) - gammaln(m) - gammaln(N - m)) / np.log(10)

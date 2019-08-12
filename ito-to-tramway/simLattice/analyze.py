"""
This file contains scripts that load, combine and analyze the trajectories generated with 'simulate.py'

    #isort:skip_file
"""

import glob
import os
import sys

import pandas as pd
from tqdm import tqdm, trange
import numpy as np

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from reinit_folder import reinit_folder


data_folder = r"D:\calculated_data\lattice_diffusivity_obstacles_N=1000"
snr_label = 'snr'
sigma = 0.00    # um
extension = '.dat'
kwargs = {
    'method': 'hexagon', 'avg_distance': 0.5, 'min_location_count': 0,
    'reset_origin': False,
    'sep': '\t',
    'xlims': [0, 10],
    'ylims': [0, 10],
    'clip_D': 0.975,
    'clip_grad': 0.9,
    'clip_alpha': 0.9,
}


all_trajectories_file = 'all_trajectories.dat'
results_folder = "bayes_factors"

max_N = 1000


def analyze(recalculate=False):
    basename = os.path.basename(os.path.normpath(data_folder))
    all_trajectories_filepath = basename + "_" + all_trajectories_file

    if (not recalculate) and os.path.exists(all_trajectories_filepath):
        df_all = pd.read_csv(all_trajectories_filepath, sep='\t')
        # print('A1', ~recalculate)
    else:
        data_files_list = [f for f in glob.iglob(
            os.path.join(data_folder, r"*" + extension), recursive=False)]

        file = data_files_list[0]
        df_all = pd.read_csv(file, sep='\t')

        for i in trange(1, np.min([max_N, len(data_files_list)]), desc='Loading data files'):
            # for file in data_files_list[1:]:  # , desc='Loading data files'):
            file = data_files_list[i]
            df = pd.read_csv(file, sep='\t')
            df_all = df_all.append(df, ignore_index=True)
        df_all.to_csv(all_trajectories_filepath, sep='\t')  # ,  na_rep='nan')

    # folder, _ = os.path.split(file_list[0])
    results_path = results_folder
    reinit_folder(results_path)

    # Process
    from calculate import calculate
    calculate(csv_file=all_trajectories_filepath, results_folder=results_path, bl_produce_maps=True,
              snr_label=snr_label, sigma=sigma, recalculate=recalculate, ticks=True, txt_extension=extension, pdf=True, png=True, **kwargs)

    # Plot all trajectories with lattice


# if __name__ == '__main__':
#     analyze(recalculate=False)

"""
A small script to recut the vlp file into regions while avoiding boundary artifacts.
Calculates relative displacements first and cuts then

compatibility - recalculate individual trajectories number for compatibility with old TRamWAy code
"""

import numpy as np
import pandas as pd

# input_file = r"./input/neuro-receptors-pre-select/old/EXP01-01-BetaS403A-PRE-PP2-15ms_cluster_112.trxyt"
# output_file = r"./input/neuro-receptors-pre-select/EXP01-01-BetaS403A-PRE-PP2-15ms_cluster_112_cut.txt"

input_file = r"./input/neuro-receptors-post-select/old/EXP01-01-BetaS403A-POST-PP2-15ms_cluster_112.trxyt"
output_file = r"./input/neuro-receptors-post-select/EXP01-01-BetaS403A-POST-PP2-15ms_cluster_112_cut.txt"

# half_width = 1
# interval = np.array([-1, 1]) * half_width
x_lims = [9.55, 11.55]
y_lims = [25.45, 26.8]
compatibility = True

original_data = pd.read_table(input_file, names=['n', 'x', 'y', 't'])
# print(original_data)

# calculate displacements
same_particle = original_data['n'].shift(-1) == original_data['n']
dx = original_data['x'].shift(-1) - original_data['x']
dy = original_data['y'].shift(-1) - original_data['y']
dt = original_data['t'].shift(-1) - original_data['t']
dx[~same_particle] = np.nan
dy[~same_particle] = np.nan
dt[~same_particle] = np.nan
original_data['dx'] = dx
original_data['dy'] = dy
original_data['dt'] = dt

# cut
cut_data = original_data[(original_data.x >= x_lims[0])
                         & (original_data.x <= x_lims[1])
                         & (original_data.y >= y_lims[0])
                         & (original_data.y <= y_lims[1])
                         & ~np.isnan(original_data.dx)
                         ]

# recalculate n for backwards compatibility
if compatibility:
    cut_data = cut_data.reset_index(drop=True)
    n = (cut_data['n'] != cut_data['n'].shift(1)) * 1
    cut_data.n = np.cumsum(n)
    print(cut_data.t)
    cut_data.t = (cut_data.t[1] - cut_data.t[0]) * cut_data.index


# save
cut_data.to_csv(output_file, sep='\t', index=False, header=True)

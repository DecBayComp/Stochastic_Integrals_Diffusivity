"""
A small script to recut the vlp file into regions while avoiding boundary artifacts.
Calculates relative displacements first and cuts then

compatibility - recalculate individual trajectories number for compatibility with old TRamWAy code
"""

import numpy as np
import pandas as pd

input_file = r"./input/vlp_full/VLP_WT_2_full.txt"
output_file = r"./input/vlp_2_2/VLP_WT_2_2_new_cut.txt"
x_lims = [6.240095, 6.240095 + 2]
y_lims = [12.940018, 12.940018 + 2]
compatibility = True

original_data = pd.read_table(input_file, names=['n', 'x', 'y', 't'])
print(original_data)

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

    cut_data.t = (cut_data.t[1] - cut_data.t[0]) * cut_data.index


# save
cut_data.to_csv(output_file, sep='\t', index=False, header=False)

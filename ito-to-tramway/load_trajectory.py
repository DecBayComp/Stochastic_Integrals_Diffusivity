"""
Load trajectories in standard or optical tweezers format.
If .rwa file present, just reload
"""

import os

import numpy as np
import pandas as pd

from tramway.core.analyses.lazy import Analyses
from tramway.helper import load_rwa, load_xyt, save_rwa


def load_trajectory(file, reload=False, reset_origin=True):
    full_name, extension = os.path.splitext(file)
    txt_file = full_name + '.txt'
    rwa_file = full_name + '.rwa'

    if reload or not os.path.isfile(rwa_file):
        # load the trajectories
        if 'Brownien' in txt_file:
            # time step (please check)
            dt = 1.0 / 65536
            # read the table
            xyt = load_xyt(txt_file, columns=['x', 'y', 'H'], reset_origin=reset_origin)
            # add missing columns
            xyt['n'] = np.ones(xyt.shape[0])
            xyt['t'] = np.arange(dt, (xyt.shape[0] + 1) * dt, dt)
        elif 'VLP' in txt_file:
            xyt = load_xyt(txt_file, columns=['n', 'x', 'y', 't',
                                              'dx', 'dy', 'dt'], reset_origin=reset_origin)
        else:
            xyt = load_xyt(txt_file, reset_origin=reset_origin)

        analysis_tree = Analyses(xyt)
        save_rwa(rwa_file, analysis_tree, force=True)
    else:
        analysis_tree = load_rwa(rwa_file)
        # print('Found in .rwa file: ')
        # print(analysis_tree)

    return analysis_tree

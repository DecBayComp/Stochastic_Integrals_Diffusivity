"""
Load trajectories in standard or optical tweezers format.
If .rwa file present, just reload
"""

import os

import numpy as np
import pandas as pd

from tramway.core.analyses.lazy import Analyses
from tramway.helper import load_rwa, load_xyt, save_rwa


def load_trajectory(file, reload=False, reset_origin=True, extension='.txt', sep='\t', dt=None, **kwargs):
    full_name, _ = os.path.splitext(file)
    txt_file = full_name + extension
    rwa_file = full_name + '.rwa'
    # print('A', extension)

    # print('TXT', txt_file, reload)

    if reload or not os.path.isfile(rwa_file):
        # load the trajectories
        if 'Brownien' in txt_file:
            # time step (please check)
            dt = 1.0 / 65536
            # read the table
            xyt = load_xyt(txt_file, columns=['x', 'y', 'H'], reset_origin=reset_origin, **kwargs)
            # add missing columns
            xyt['n'] = np.ones(xyt.shape[0])
            xyt['t'] = np.arange(dt, (xyt.shape[0] + 1) * dt, dt)
            xyt.drop(columns='H', inplace=True)
        elif 'VLP' in txt_file:
            xyt = load_xyt(txt_file, columns=['n', 'x', 'y', 't',
                                              'dx', 'dy', 'dt'], reset_origin=reset_origin, sep=sep)
        elif '.trxyt' in txt_file:
            # print('Good')
            print(extension, txt_file)
            xyt = load_xyt(txt_file, columns=['n', 'x', 'y', 't'],
                           reset_origin=reset_origin, sep=sep)

        elif 'EXP' in txt_file:
            # print('Good')
            print(extension, txt_file)
            xyt = load_xyt(txt_file, columns=['n', 'x', 'y', 't',
                                              'dx', 'dy', 'dt'], reset_origin=reset_origin, sep=sep)

        elif 'sim_data' in txt_file:
            print(extension, txt_file)
            xyt = load_xyt(txt_file, columns='x dx y dy'.split(),
                           reset_origin=reset_origin, sep=sep)
            xyt['n'] = np.ones(xyt.shape[0])
            xyt['t'] = np.arange(dt, (xyt.shape[0] + 1) * dt, dt)

        #     # print('Good')
        #     print(extension, txt_file)
        #     xyt = pd.read_csv(txt_file, sep=';', skiprows=1,
        #                       names=['x', 'dx', 'y', 'dy'])
        #     xyt.drop(xyt.tail(1).index, inplace=True)
        #     xyt['n'] = np.zeros(xyt.shape[0])
        #     xyt['t'] = np.arange(0., dt * xyt.shape[0], dt)
        #     xyt['dt'] = dt

        # elif '.trxyt' in txt_file:
        #     # print('Good')
        #     print(extension, txt_file)
        #     xyt = load_xyt(txt_file, columns=['n', 'x', 'y', 't'], reset_origin=reset_origin)

            # print(type(xyt))
            # print(xyt)
        else:
            xyt = load_xyt(txt_file, reset_origin=reset_origin, sep=sep)

        # Restore the dx, dy, dt columns if absent
        # print(xyt)
        # print(xyt.keys())
        # print(xyt.columns.values)
        keys = (key for key in ['x', 'y', 't'] if 'd' + key not in xyt.columns.values)
        ns = xyt.n.unique()
        for key in keys:
            xyt['d' + key] = np.nan
            for n in ns:
                inds = xyt.n == n
                xyt.loc[inds, 'd' + key] = xyt.loc[inds,
                                                   key].shift(periods=-1) - xyt.loc[inds, key]

        # print('C', xyt)

        analysis_tree = Analyses(xyt)
        save_rwa(rwa_file, analysis_tree, force=True)
    else:
        analysis_tree = load_rwa(rwa_file)
        # print('Found in .rwa file: ')
        # print(analysis_tree)

    return analysis_tree

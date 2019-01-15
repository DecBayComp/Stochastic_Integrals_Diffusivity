#!/usr/bin/env python

"""
This function may be outdated
"""

import os.path

import numpy as np
import pandas as pd

from constants import dt
from Sashas import Sashas
from snr import infer_snr
from tramway.core import *
from tramway.helper import *
from tramway.helper import Analyses, infer
from tramway.inference import *


def tesselate_and_infer(csv_file, bl_dr=True, localization_error=1e-3, load=True):
    rwa_file, _ = os.path.splitext(csv_file)
    rwa_file = '{}.rwa'.format(rwa_file)

    if not os.path.exists(rwa_file) or not load:
        data = pd.read_csv(csv_file, sep=';', skiprows=1,
                           names=['x', 'dx', 'y', 'dy'])
        data.drop(data.tail(1).index, inplace=True)
        data['n'] = np.zeros(data.shape[0])
        data['t'] = np.arange(0., dt * data.shape[0], dt)
        data['dt'] = dt
        # print(data)

        analysis_tree = Analyses(data)
        tessellate(analysis_tree, 'gwr', output_file=rwa_file, label='gwr',
                   verbose=True, force=True)  # , ref_distance = 0.05)
        # cell_plot(rwa_file, output_file = traj+'_mesh.png')
    else:
        print("Warning: trajectories not reloaded! Meshes not recalculated.")

    print("\nMesh constructed. Starting inference\n")
    # rwa_file = '.\\results\\{}.rwa'.format(traj)
    # a = load_rwa(rwa_file)
    if 'dx' in analysis_tree.data:
        cell_types = Sashas
    else:
        cell_types = Translocations
    # print(analysis_tree.data)
    snr_label = 'snr'
    infer(analysis_tree, 'snr', input_label='gwr', output_label=snr_label,
          max_iter=50, sigma=0, new_cell=cell_types)

    # d = distributed(analysis_tree['gwr'].data, new_cell=cell_types)
    # # , diffusivity_prior = 0.001, min_diffusivity = 0.001)
    # x = d.run(infer_snr, max_iter=50, localization_error=localization_error)
    # a['gwr'].add(Maps(x, mode='snr'), label='snr')
    save_rwa(rwa_file, analysis_tree, verbose=True, force=True)

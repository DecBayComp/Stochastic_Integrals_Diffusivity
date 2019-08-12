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


def tesselate_and_infer(csv_file, bl_dr=True, sigma2=1e-3, load=True, min_diffusivity=1e-8, **kwargs):
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
                   verbose=True, force=True, **kwargs)  # , ref_distance = 0.05)
    else:
        print("Warning: trajectories not reloaded! Meshes not recalculated.")
        analysis_tree = load_rwa(rwa_file)

    print("\nMesh constructed. Starting inference\n")

    snr_label = 'snr'
    infer(analysis_tree, 'd.conj_prior', input_label='gwr', output_label=snr_label,
          sigma2=sigma2, **kwargs)

    save_rwa(rwa_file, analysis_tree, verbose=True, force=True)

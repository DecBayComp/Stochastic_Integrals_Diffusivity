#!/usr/bin/env python

from tramway.core import *
from tramway.helper import *
from tramway.inference import *
import numpy as np
import pandas as pd
import os.path
from snr import infer_snr
from constants import dt



def tesselate_and_infer(csv_file, bl_dr = True):
    rwa_file, _ = os.path.splitext(csv_file)
    rwa_file = '{}.rwa'.format(rwa_file)

    if not os.path.exists(rwa_file):
        data = pd.read_csv(csv_file, sep=';', skiprows=1, names=['x', 'dx', 'y', 'dy'])
        data['n'] = np.zeros(data.shape[0])
        data['t'] = np.arange(0., dt * data.shape[0], dt)

        tessellate(data, 'gwr', output_file=rwa_file, label='gwr', verbose=True, force = True) #, ref_distance = 0.05)
        # cell_plot(rwa_file, output_file = traj+'_mesh.png')
    else:
        print("Warning: trajectories not reloaded! Meshes not recalculated.")

    class Sashas(Translocations):
        def _extract_space(self):
            return np.asarray(self.origins[['dx', 'dy']])
        @property
        def space_cols(self):
            return ['x', 'y']
        @space_cols.setter
        def space_cols(self, cs):
            pass

    print("\nMesh constructed. Starting inference\n")
    # rwa_file = '.\\results\\{}.rwa'.format(traj)
    a = load_rwa(rwa_file)
    if bl_dr:
        cell_types = Sashas
    else:
        cell_types = Translocations
    d = distributed(a['gwr'].data, new_cell=cell_types)
    x = d.run(infer_snr, max_iter=50, localization_error = 1e-3) #, diffusivity_prior = 0.001, min_diffusivity = 0.001)
    a['gwr'].add(Maps(x, mode='snr'), label='snr')
    save_rwa(rwa_file, a, verbose = True, force = True)


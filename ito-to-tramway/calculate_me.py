#!/usr/bin/env python

from tramway.core import *
from tramway.helper import *
from tramway.inference import *
import numpy as np
import pandas as pd
import os.path
from snr import infer_snr

dt = .04

traj_file = "./data/000000002.csv"
traj, _ = os.path.splitext(traj_file)
traj = traj.split('_')[-1]
rwa_file = '{}.rwa'.format(traj)

if not os.path.exists(rwa_file):
	data = pd.read_csv(traj_file, sep=';', skiprows=1, names=['x', 'dx', 'y', 'dy'])
	data['n'] = np.zeros(data.shape[0])
	data['t'] = np.arange(0., dt * data.shape[0], dt)

	tessellate(data, 'gwr', output_file=rwa_file, label='gwr', verbose=True)
	cell_plot(rwa_file, output_file=traj+'.png')

# from tramway.inference import *
# setup, module = plugins['snr']
# snr = getattr(module, setup['infer'])

class Sashas(Translocations):
	def _extract_space(self):
		return np.asarray(self.origins[['dx', 'dy']])
	@property
	def space_cols(self):
		return ['x', 'y']
	@space_cols.setter
	def space_cols(self, cs):
		pass

# rwa_file = '.\\results\\{}.rwa'.format(traj)
a = load_rwa(rwa_file)
d = distributed(a['gwr'].data, new_cell=Sashas)
x = d.run(infer_snr, max_iter=50, localization_error = 0.0)
a['gwr'].add(Maps(x, mode='snr'), label='snr')
save_rwa(rwa_file, a)


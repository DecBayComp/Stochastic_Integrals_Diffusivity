

from constants import optical_traps_data_folder, optical_data_sets, optical_traps_dt
import numpy as np
import os
import pandas as pd


def load_trajectory(filename):

    # Init
    input_folder = optical_traps_data_folder

    txt_file = filename + '.txt'
    rwa_file = filename + '.rwa'

    txt_fullpath = os.path.join(input_folder, txt_file)
    rwa_fullpath = os.path.join(input_folder, rwa_file)

    trajectory = pd.read_csv(txt_fullpath, sep='\t', names=['x', 'y', 'smth'])
    del trajectory['smth']

    # Reset mean x and y to 0
    trajectory['x'] -= trajectory['x'].mean()
    trajectory['y'] -= trajectory['y'].mean()

    # Assign track number and time
    trajectory['n'] = np.zeros(trajectory.shape[0])
    trajectory['t'] = np.arange(0., optical_traps_dt * trajectory.shape[0], optical_traps_dt)

    return [trajectory, rwa_fullpath]

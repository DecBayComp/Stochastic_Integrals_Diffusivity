"""
Parameters for batch trajectories analysis with TRamWAy and Bayes factor plots
"""

import os  # the right way to join paths
import platform

import numpy as np

version = 20180416

# Time step for inference
dt = 0.04  # s
L = 1  # um

# Simulation parameters for comparison
D_0 = 0.01  # um^2/s
D_ratio = 2.0
zeta_t_y_over_zeta_sp_abs = 6.25
k = 2.0  # um^{-1}, D'/D_0
# ksi = -10.0
abs_tol = 1.0e-8
bl_produce_maps = False
# remember to change the value for the analyzed data. For theory use a value much smaller than the jump
sigma = abs_tol


# if sys ==

# data_folder_lin = '/mnt/d/calculated_data/sim_performance_2D_no_perp'
# data_folder_win = r'd:\calculated_data\sim_performance_2D_no_perp'
# data_folder_lin = '/mnt/d/calculated_data/sim_performance_2D_with_perp'
# data_folder_win = r'd:\calculated_data\sim_performance_2D_with_perp'

# Neuro-receptors data
data_folder_lin = '/mnt/d/Google Drive/git/Stochastic_Integrals_Diffusivity/ito - to - tramway/input/neuro - receptors'
data_folder_win = r'd:\Google Drive\git\Stochastic_Integrals_Diffusivity\ito - to - tramway\input\neuro - receptors'

args_file = "./arguments.dat"
args_lock = "./arguments.lock"

optical_traps_data_folder_win = r'D:\Experimental_Data\optical_tweezers'
optical_traps_data_folder_lin = '/mnt/d/Experimental_Data/optical_tweezers'


optical_data_sets = ['P1', 'P4', 'P5']
optical_power_mW = [500, 251, 138]
optical_traps_points_per_bin = 400
optical_traps_dt = 1.0 / 65536


def folders():
    sys = platform.system()
    if sys == "Windows":
        data_folder = data_folder_win
        optical_traps_data_folder = optical_traps_data_folder_win
    elif sys == "Linux" or sys == "Darwin":
        data_folder = data_folder_lin
        optical_traps_data_folder = optical_traps_data_folder_lin
    else:
        print("Error: unknown OS (%i)" % sys)
    output_folder = os.path.join(data_folder, "dat")

    return data_folder, output_folder


# Identify the optical data traps folder
sys = platform.system()
if sys == "Windows":
    optical_traps_data_folder = optical_traps_data_folder_win
elif sys == "Linux" or sys == "Darwin":
    optical_traps_data_folder = optical_traps_data_folder_lin
else:
    print("Error: unknown OS (%i)" % sys)


CSV_DELIMITER = ';'

# Batch parameters
trials = 1  # 1000
sleep_time = 0.2
logs_folder = "./logs/"
lock_timeout = 300

jobs_count_tars = 1  # 132+12
jobs_count_t_bayes = 132
jobs_count_onsager = 10  # 10
manager_script = "job_manager.py"
DETACHED_PROCESS = 0x00000008

# # >> round-well calculations <<
# data_folder = './data/round_well/'
# output_folder = './results/round_well/'
# args_file = "./arguments_2.dat"
# args_lock = "./arguments_2.lock"

# color sequence
alpha = 0.25
green = [0.4353,    0.5804,         0]
red = np.asarray([191,    32,    37]) / 255.0
yellow = np.asarray([252.91, 191.76, 16.47]) / 255.0

col_yes = np.asarray([194,    216, 52]) / 255.0
col_no = np.asarray([141,  24, 26]) / 255.0
col_idk = np.asarray([255,   255,    255]) / 255.0
colors = np.transpose(np.stack([col_no, col_idk, col_yes], axis=1))

MACHINE_PRECISION = 1e-16
pagewidth_in = 6.85

color_sequence = [[0.1008,    0.4407,    0.7238],
                  [0.4353,    0.5804,         0],
                  [0.9498,    0.4075,    0.1317],
                  [0.4364,    0.2238,    0.5872],
                  [0.5860,    0.4228,    0.2649]]

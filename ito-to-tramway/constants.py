"""
Parameters for batch trajectories analysis with TRamWAy and Bayes factor plots
"""

import numpy as np
import os   # the right way to join paths
import platform

version = 20180416

# Time step for inference
dt = 0.04  # s
L = 1  # um

# Simulation parameters for comparison
D_0 = 0.01  # um^2/s
D_ratio = 2.0
k = 2.0  # um^{-1}, D'/D_0
# ksi = -10.0
abs_tol = 1.0e-8
bl_produce_maps = False
# remember to change the value for the analyzed data. For theory use a value much smaller than the jump
localization_error = abs_tol


# if sys ==

# data_folder_lin = '/mnt/d/calculated_data/sim_performance_2D_no_perp'
# data_folder_win = r'd:\calculated_data\sim_performance_2D_no_perp'
data_folder_lin = '/mnt/d/calculated_data/sim_performance_2D_with_perp'
data_folder_win = r'd:\calculated_data\sim_performance_2D_with_perp'

args_file = "./arguments.dat"
args_lock = "./arguments.lock"


def folders():
    sys = platform.system()
    if sys == "Windows":
        data_folder = data_folder_win
    elif sys == "Linux" or sys == "Darwin":
        data_folder = data_folder_lin
    else:
        print("Error: unknown OS (%i)" % sys)
    output_folder = os.path.join(data_folder, "dat")

    return data_folder, output_folder


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
green = [0.4353,    0.5804,         0]
red = np.asarray([191,    32,    37]) / 255.0
yellow = np.asarray([252.91, 191.76, 16.47]) / 255.0

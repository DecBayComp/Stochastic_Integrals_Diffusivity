"""
Parameters for batch trajectories analysis with TRamWAy and Bayes factor plots
"""

version = 20180416

# Time step for inference
dt = 0.04  # s

# Simulation parameters for comparison
D_0 = 0.01  # um^2/s
D_ratio = 2.0
k = 2.0  # um^{-1}
# ksi = -10.0
abs_tol = 1.0e-8
bl_produce_maps = True

# >> main validation with the saw-tooth profile <<
data_folder = '/mnt/d/calculated_data/bayes_test_tramway_2D/'
output_folder = './results/'
args_file = "./arguments.dat"
args_lock = "./arguments.lock"

CSV_DELIMITER = ';'

# Batch parameters
trials = 1  # 1000
sleep_time = 0.2
logs_folder = "./logs/"
lock_timeout = 300

jobs_count_tars = 1  # 132+12
jobs_count_t_bayes = 132
jobs_count_onsager = 12
manager_script = "job_manager.py"
DETACHED_PROCESS = 0x00000008

# # >> round-well calculations <<
# data_folder = './data/round_well/'
# output_folder = './results/round_well/'
# args_file = "./arguments_2.dat"
# args_lock = "./arguments_2.lock"

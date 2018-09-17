"""
This script collects .rwa files in subfolders on the network and creates an arguments file with their names.
The goal is to analyze J-B's simulated trajectories
"""

try:
    bl_has_run
except:
    %load_ext autoreload
    %autoreload 2
    bl_has_run = True

from calculate import calculate
import glob
import matplotlib
import os

# % Constants
# dt = 0.05
# snr_label = 'snr(mu=0)'
# snr_label = 'snr'
results_folder = "bayes_factors"


# root_path = r"\\157.99.40.171\@Dbc\LAB_shared_stuff\Francois_Laurent\tests_tramway\numerical_trajectories_no_box\snr"
# root_path = r"\\157.99.40.171\@Dbc\LAB_shared_stuff\Francois_Laurent\tests_tramway\numerical_trajectories"

# # Miscellaneous experiments
# root_path = r"\\157.99.40.171\@Dbc\LAB_shared_stuff\Francois_Laurent\tests_tramway\misc_experiments"
# dt = 0.05  # s
# snr_label = 'snr'
# localization_error = 0.03
# bl_recursive = False

# # optical tweezers
# root_path = r"\\157.99.40.171\@Dbc\LAB_shared_stuff\Francois_Laurent\tests_tramway\optical_tweezers"
# dt = 1 / 65636  # s
# snr_label = 'snr(mu=0)'
# localization_error = 0.01

# # Full viral capside
# root_path = r"\\157.99.40.171\@Dbc\LAB_shared_stuff\Francois_Laurent\tests_tramway\misc_experiments\full"
# dt = 0.02  # s
# snr_label = 'snr'
# localization_error = 0.03
# bl_recursive = False

# VLP region
root_path = r"\\157.99.40.171\@Dbc\LAB_shared_stuff\Francois_Laurent\tests_tramway\misc_experiments\vlp_2_2"
dt = 0.02  # s
snr_label = 'snr'
localization_error = 0.03
bl_recursive = False

# root_path = r"D:\\git\Stochastic_Integrals_Diffusivity\ito-to-tramway"


# %% Parse all subdirectories to extract .rwa files and store
# print(os.path.exists(root_path))
# print(root_path + r"\**\*.rwa")

if bl_recursive:
    file_list = [f for f in glob.iglob(root_path + r"\**\*.rwa", recursive=True)]
else:
    file_list = [f for f in glob.iglob(root_path + r"\*.rwa", recursive=False)]

# print("Found files: \n", file_list)
# print(file_list)

for i in range(len(file_list)):
    file = file_list[i]
    folder, _ = os.path.split(file)
    folder = folder + "\\"
    print("Processing file: %s" % file)
    # calculate(file, folder, True, dt, snr_label)
    results_path = os.path.join(folder, results_folder)
    calculate(csv_file=file, results_folder=results_path, bl_produce_maps=True,
              dt=dt, snr_label=snr_label, localization_error=localization_error)

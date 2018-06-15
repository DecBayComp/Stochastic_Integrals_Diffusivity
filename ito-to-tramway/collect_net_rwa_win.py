"""
This script collects .rwa files in subfolders on the network and creates an arguments file with their names.
The goal is to analyze J-B's simulated trajectories
"""

%load_ext autoreload
%autoreload 2

from calculate import calculate
import glob
import os

# % Constants
dt = 0.05
snr_label = 'snr(mu=0)'
snr_label = 'snr'
results_folder = "bayes_factors"
localization_error = 0.01


# root_path = r"\\157.99.40.171\@Dbc\LAB_shared_stuff\Francois_Laurent\tests_tramway\numerical_trajectories_no_box\snr"
# root_path = r"\\157.99.40.171\@Dbc\LAB_shared_stuff\Francois_Laurent\tests_tramway\numerical_trajectories"
root_path = r"\\157.99.40.171\@Dbc\LAB_shared_stuff\Francois_Laurent\tests_tramway\misc_experiments"
# root_path = r"D:\\git\Stochastic_Integrals_Diffusivity\ito-to-tramway"


# %% Parse all subdirectories to extract .rwa files and store
# print(os.path.exists(root_path))
# print(root_path + r"\**\*.rwa")
file_list = [f for f in glob.iglob(root_path + r"\**\*.rwa", recursive=True)]
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

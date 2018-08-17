"""
Calculate and plot the Bayes factors in optical traps provided by Max. The idea is to keep the number of jumps per bin the same, but change the power of the optical trap. Due to this, there will be different level of evidence for the active forces in the different cases, although in all cases there is confinement. In the weakest cases, the confienement might be due to another stochastic integral interpretation.
"""

try:
    has_run
except NameError:
    %matplotlib
    %load_ext autoreload
    %autoreload 2

    has_run = 1
else:
    print("Graphic interface NOT re-initialized")

from calculate import calculate
from constants import optical_traps_data_folder, optical_data_sets, optical_traps_points_per_bin, MACHINE_PRECISION, optical_traps_dt, localization_error
import matplotlib
import matplotlib.pyplot as plt
import os
import pandas as pd
from snr import infer_snr
from tramway.helper import *

# def plot_optical_traps():

# %% Constants
font_size = 8
alpha = 1.0
linewidth = 0.1
dpi = 100

# %% Load and create a uniform jump number mesh for some optical trap data files
input_folder = optical_traps_data_folder
file_ind = 1
txt_file = 'Brownien_%s.txt' % (optical_data_sets[file_ind])
txt_fullpath = os.path.join(input_folder, txt_file)

rwa_file = 'Brownien_%s.rwa' % (optical_data_sets[file_ind])
rwa_fullpath = os.path.join(input_folder, rwa_file)

trajectory = pd.read_csv(txt_fullpath, sep='\t', names=['x', 'y', 'smth'])
del trajectory['smth']
# Assign track number and time
trajectory['n'] = np.zeros(trajectory.shape[0])
trajectory['t'] = np.arange(0., optical_traps_dt * trajectory.shape[0], optical_traps_dt)

# %% Tessellate with a fixed number of jumps per bin
# mesh_ = tessellate(trajectory, 'hexagon', verbose = True, avg_location_count=count, ref_distance=MACHINE_PRECISION, force=True)
mesh = tessellate(trajectory, 'kmeans', label='equal', verbose = True, avg_location_count=optical_traps_points_per_bin, ref_distance=MACHINE_PRECISION, prune = 0, force=True)


fig = plt.figure(1)
plt.clf()
cell_plot(mesh, num = 1, figsize = None)

# Initialize the tree
rwa_data = Analyses(trajectory)
# Add a mesh
rwa_data.add(mesh, label='equal')
equal_tree = rwa_data['equal']

# Change the cell representation
equal_tree_distr = distributed(equal_tree.data)
ns = np.asarray([len(equal_tree_distr[i]) for i in equal_tree_distr])
print(ns)
 
# %% Perform inference and calculate the Bayes factors in each bin
# -> May need to check the localization error value
print(rwa_data)
snr_result = equal_tree_distr.run(infer_snr, max_iter = 50, localization_error = localization_error)
equal_tree.add(Maps(snr_result, mode = 'snr'), label = 'snr')
print(rwa_data)
save_rwa(rwa_fullpath, rwa_data, verbose=4, force=True)

# %% Calculate the Bayes factors
output_folder = ".\output"
calculate(rwa_fullpath, output_folder, bl_produce_maps=False, dt=optical_traps_dt,
          snr_label='snr', localization_error=localization_error)

# Reload to get the calculated Bayes factors
rwa_data = load_rwa(rwa_fullpath)
cells = rwa_data['equal'].data
n = rwa_data['equal']['snr'].data['n']
log10_Bs = rwa_data['equal']['snr']['bayes_factors'].data.log10_B

# Plot a map of the Bayes factors
my_map = pd.DataFrame(np.asarray(log10_Bs), index=n.index, columns=['$\log_{10}(B)$'])
# fig = plt.figure(2)
matplotlib.rcParams.update({'font.size': font_size})
map_plot(my_map, cells=cells, show=False, clip=False, colormap='inferno', alpha=alpha, linewidth=linewidth, dpi=dpi, aspect='equal') #, labelsize=font_size)
# fig.show()
         
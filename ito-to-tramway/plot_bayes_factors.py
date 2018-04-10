#!/usr/bin/env python

from tramway.core   import *
from tramway.helper import *
import os.path
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt

from tramway.inference.calculate_bayes_factors import calculate_bayes_factors


# all the following examples are 2D
precomputed_meshes = [
	# standard example where nothing special happens
	'glycine_receptor',
	# lipid raft example from lipid_platform/wild/folder_2008-07-04-lamella1_NP-PTE_20mW-50ms_Series14_T27/2103.txt
	'lipid_2103',
	# lipid raft example from lipid_platform/wild/folder_2012-02-29_CS1_20mW_50ms_Series07/3952.txt
	'lipid_3952',
	# lipid raft example from lipid_platform/transferin/folder_2014-12-17_30mW_50ms_SilvanNP_CS1_Series06/4305.txt
	'lipid_4305',
	# VLP example from VLP/WT/2/trajectories_2.txt
	'VLP_WT_2_2',
	]

# my advice: start with a single mesh
precomputed_meshes = [ '000000002' ]

# lipid_* and VLP_* files contain several kmeans and gwr meshes of varying spatial resolution;
# label `mesh` is either 'kmeans' or 'gwr' followed by a number that approximates the average 
# location count per cell (20, 40 or 60)


snr_label = 'snr'
Vpi_name = 'V_prior'

# a few convenience functions
def to_array(a):
	a = a.values
	return a if a.shape[1:] else a[:,np.newaxis]
sum_cells = lambda a: np.sum(a, axis=0, keepdims=True)
sum_dims  = lambda a: np.sum(a, axis=1, keepdims=True)
vec_to_2D = lambda a: (np.asarray(a) @ (np.asarray([[1, 1]])))

# loop over the selected examples
for example in precomputed_meshes:

	# ensure the .rwa file exists
	# (will skip to the next step for precomputed meshes/maps)
	if not os.path.exists(example+'.rwa'):
		if not os.path.isfile(example+'.txt'):
			raise OSError("cannot find file '{}.[txt|rwa]'".format(example))
		# generate a default mesh
		tessellate(example+'.txt', 'gwr', strict_min_location_count=10, label='gwr')

	# load the .rwa file
	analysis_tree = load_rwa(example+'.rwa')

	# loop over the available meshes
	anything_new = False
	for mesh in analysis_tree: # `mesh` is a label (a key in dict-like `analysis_tree`)
		print('{} - {}'.format(example, mesh))

		# infer snr-related maps with the 'snr' plugin
		# (will skip to the next step for precomputed maps)
		if snr_label not in analysis_tree[mesh]:
			# `infer` adds a subtree identified by `snr_label`
			infer(analysis_tree, 'snr', input_label=mesh, output_label=snr_label,
				max_iter=50)
			anything_new = True
			save_rwa(example+'.rwa', analysis_tree, force=True)

		snr = analysis_tree[mesh][snr_label].data

		# individualize the various intermediate maps, namely:
		# - jump count `n`,
		# - variance `V_prior` of all the jumps but those in the cell,
		# - total snr `zeta_total`,
		# - spurious snr `zeta_spurious`.

		# print(type(snr))
		n = snr['n']
		D = snr['diffusivity']
		print(np.asarray(D))
		zeta_total = snr['zeta_total']
		zeta_spurious = snr['zeta_spurious']
		if True: # Vpi_name not in snr.variables:
			# `V_prior` is not computed directly by the 'snr' plugin
			# because the plugin's main routine may be independently applied to 
			# subsets of cells instead of all the cells at a time
			dr, dr2   = to_array(snr['dr']), to_array(snr['dr2'])
			n_prior   = np.sum(to_array(n))  - to_array(n)
			dr_prior  = sum_cells(dr) - dr
			dr2_prior = sum_cells(dr2) - dr2
			
			
			# calculate biased varaince in current bin
			dr_mean = dr / vec_to_2D(n)
			dr_mean2 = sum_dims(dr_mean ** 2)
			dr2_mean = sum_dims(dr2) / n
			Vs = np.asarray(dr2_mean - dr_mean2)

			# calculate prior variance
			dr_prior_mean = dr_prior / n_prior
			dr_prior_mean2 = sum_dims(dr_prior_mean ** 2)
			dr2_prior_mean = sum_dims(dr2_prior) / n_prior
			Vs_prior = np.asarray(dr2_prior_mean - dr_prior_mean2)


			# print(Vs_prior / Vs)

			# V = 

			# V_prior   = sum_dims(dr2_prior - dr_prior * dr_prior / n_prior) / (n_prior - 1)
			# add `V_prior` to the analysis tree before the latter is saved
			# snr.maps  = snr.maps.join(pd.DataFrame(
			# 	data    = V_prior,
			# 	index   = n.index,
			# 	columns = [Vpi_name],
			# 	))
			# anything_new = True
		# V_prior = snr[Vpi_name]

		# TODO: generate final maps, say the `my_map` map
		
		# Prepare numpy arrays for Bayes factor calculation
		zeta_ts = np.asarray(zeta_total)
		zeta_sps = np.asarray(zeta_spurious)
		ns = np.asarray(n)
		D = np.asarray(D)
		# Vs = 

		Bs, bl_forces = calculate_bayes_factors(zeta_ts = zeta_ts, zeta_sps = zeta_sps, ns = ns, Vs = Vs, Vs_pi = Vs_prior)
		# print (bl_forces)

		

		# Check the inferred gradient compared to the simulated gradient
		dt = 0.04 # s
		D_0 = 0.01 # um^2/s
		k = 2.0 # um^{-1}
		ksi = -1.0
		abs_tol = 1.0e-8
		
		# Set non-calculated values to nan
		zeta_sps[np.abs(zeta_sps) < abs_tol] = np.nan
		# print(zeta_sps)
		gdt_inf = zeta_sps * vec_to_2D(np.sqrt(Vs))
		mean_gdt_inf = np.nanmedian(gdt_inf, axis = 0)
		gdt_sim = np.asarray([D_0 * k, 0]) * dt

		# Check total force
		alpha_dt_inf = zeta_ts * vec_to_2D(np.sqrt(Vs))
		mean_alpha_dt_inf = np.median(alpha_dt_inf, axis = 0)
		alpha_dt_sim = np.asarray([ksi * D_0 * k , 0]) * dt
		
		# Simulated SNRs
		zeta_ts_exp = alpha_dt_sim / vec_to_2D(np.sqrt(Vs))
		zeta_sps_exp = gdt_sim / vec_to_2D(np.sqrt(Vs))

		# General parameters
		print ("\n\nMedian Vs:\t%.2g um^2" % (np.median(Vs)))
		print ("Median sqrt(Vs):\t%.2g um" % (np.median(np.sqrt(Vs))))
		print ("Median zeta_ts:\t%s" % (np.median(zeta_ts, axis = 0)))
		print ("Expected zeta_ts:\t%s" % (np.median(zeta_ts_exp, axis = 0)))
		print ("Median zeta_sps:\t%s" % (np.nanmedian(zeta_sps, axis = 0)))
		print ("Expected zeta_sps:\t%s" % (np.median(zeta_sps_exp, axis = 0)))
		print ("Median n:\t%i" % np.median(ns, axis = 0))

		# Inferred gradient
		print("Inferred <g dt>:\t%s um" % (mean_gdt_inf))
		print("Simulated <g dt>:\t%s um" % (gdt_sim))

		# Check total force
		print("Inferred <alpha dt>:\t%s um" % (mean_alpha_dt_inf))
		print("Simulated <alpha dt>:\t%s um" % (alpha_dt_sim))

		print("\n\n")

		# my_map = pd.DataFrame(np.log10(Bs), index = n.index, columns = ['log10(B)'])
		# my_map = pd.DataFrame(bl_forces * 1.0, index = n.index, columns = ['bool Bayes factor'])
		my_map = pd.DataFrame(D.T[0], index = n.index, columns = ['D'])
		# my_map = pd.DataFrame(alpha_dt_inf, index = n.index, columns = ['$alpha dt$ x', '$alpha dt$ y'])
		# my_map = pd.DataFrame(gdt_inf, index = n.index, columns = ['$g dt$ x', '$g dt$ y'])

		# the above line is intended to be replaced by more complex calculations;
		# it ensures the following lines do not crash because of missing variable `my_map`

		# plot the `my_map` map
		cells = analysis_tree[mesh].data # `cells` contains the mesh
		# map_plot(my_map, cells=cells, show=False)
		# ...
		# plt.show() # waits for the user to close the resulting window
		## ... or alternatively plot in a file
		map_plot(my_map, cells=cells, output_file = 'result_' + precomputed_meshes[0] + "_" + mesh + '.png', clip = False)
		


	# save the intermediate snr-related maps
	if anything_new:
		save_rwa(example+'.rwa', analysis_tree, force=True)


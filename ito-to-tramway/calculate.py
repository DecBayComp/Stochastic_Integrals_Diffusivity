#!/usr/bin/env python

"""
Calculate Bayes factors after TRamWAY basic inference if finished.
Plots figures.
"""

import csv
from tramway.core import *
from tramway.helper import *
import os.path
import numpy as np
import pandas as pd
import matplotlib
# matplotlib.use('Agg') # enable for console runs with no displays
import matplotlib.pyplot as plt

from bayes_factors.calculate_bayes_factors import calculate_bayes_factors
from constants import abs_tol, data_folder, CSV_DELIMITER

# for theoretical checks
from constants import D_0, k


def calculate(csv_file, results_folder, bl_output_map, dt, snr_label):
    # # all the following examples are 2D
    # precomputed_meshes = [
    # 	# standard example where nothing special happens
    # 	'glycine_receptor',
    # 	# lipid raft example from lipid_platform/wild/folder_2008-07-04-lamella1_NP-PTE_20mW-50ms_Series14_T27/2103.txt
    # 	'lipid_2103',
    # 	# lipid raft example from lipid_platform/wild/folder_2012-02-29_CS1_20mW_50ms_Series07/3952.txt
    # 	'lipid_3952',
    # 	# lipid raft example from lipid_platform/transferin/folder_2014-12-17_30mW_50ms_SilvanNP_CS1_Series06/4305.txt
    # 	'lipid_4305',
    # 	# VLP example from VLP/WT/2/trajectories_2.txt
    # 	'VLP_WT_2_2',
    # 	]

    # # my advice: start with a single mesh
    # precomputed_meshes = [ data_folder + filename ]

    # # lipid_* and VLP_* files contain several kmeans and gwr meshes of varying spatial resolution;
    # # label `mesh` is either 'kmeans' or 'gwr' followed by a number that approximates the average
    # # location count per cell (20, 40 or 60)

    # Replace .csv extension by .rwa
    input_file_no_ext, _ = os.path.splitext(csv_file)
    rwa_file = '{}.rwa'.format(input_file_no_ext)
    filename = os.path.splitext(os.path.basename(csv_file))[0]
    if not os.path.exists(rwa_file):
        raise RuntimeError("rwa file (%s) not found. Aborting" % rwa_file)

    # Read ksi value from .csv file if it exists
    ksi = np.nan
    try:
        with open(csv_file, 'r') as f:
            reader = csv.reader(f, delimiter=CSV_DELIMITER)
            ksi = float(next(reader)[2])
            # print(ksi)
    except:
        print(".csv file not found. Skipping reading ksi")

    # snr_label = 'snr(mu=0)'
    Vpi_name = 'V_prior'

    # a few convenience functions
    def to_array(a):
        a = a.values
        return a if a.shape[1:] else a[:, np.newaxis]

    def sum_cells(a): return np.sum(a, axis=0, keepdims=True)

    def sum_dims(a): return np.sum(a, axis=1, keepdims=True)

    def vec_to_2D(a): return (np.asarray(a) @ (np.asarray([[1, 1]])))

    # # loop over the selected examples
    # for example in precomputed_meshes:

    # # ensure the .rwa file exists
    # # (will skip to the next step for precomputed meshes/maps)
    # if not os.path.exists(example+'.rwa'):
    # 	if not os.path.isfile(example+'.txt'):
    # 		raise OSError("cannot find file '{}.[txt|rwa]'".format(example))
    # 	# generate a default mesh
    # 	tessellate(example+'.txt', 'gwr', strict_min_location_count=10, label='gwr')

    # load the .rwa file
    analysis_tree = load_rwa(rwa_file)

    # loop over the available meshes
    anything_new = False
    # `mesh` is a label (a key in dict-like `analysis_tree`)
    for mesh in analysis_tree:
        print("Using ", rwa_file)

        # infer snr-related maps with the 'snr' plugin
        # (will skip to the next step for precomputed maps)
        if snr_label not in analysis_tree[mesh]:
            # `infer` adds a subtree identified by `snr_label`
            infer(analysis_tree, 'snr', input_label=mesh, output_label=snr_label,
                  max_iter=50)
            anything_new = True
            save_rwa(rwa_file, analysis_tree, force=True)

        snr = analysis_tree[mesh][snr_label].data

        # individualize the various intermediate maps, namely:
        # - jump count `n`,
        # - variance `V_prior` of all the jumps but those in the cell,
        # - total snr `zeta_total`,
        # - spurious snr `zeta_spurious`.

        # print(type(snr))
        n = snr['n']
        D = snr['diffusivity']
        # print(np.asarray(D))
        zeta_total = snr['zeta_total']
        zeta_spurious = snr['zeta_spurious']
        if True:  # Vpi_name not in snr.variables:
            # `V_prior` is not computed directly by the 'snr' plugin
            # because the plugin's main routine may be independently applied to
            # subsets of cells instead of all the cells at a time
            dr, dr2 = to_array(snr['dr']), to_array(snr['dr2'])
            n_prior = np.sum(to_array(n)) - to_array(n)
            dr_prior = sum_cells(dr) - dr
            dr2_prior = sum_cells(dr2) - dr2

            # calculate biased varaince in current bin
            dr_mean = dr / vec_to_2D(n)
            dr_mean2 = sum_dims(dr_mean ** 2)
            dr2_mean = sum_dims(dr2) / n
            # print(np.mean(dr_mean2))
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

        # print("Vs: %s" % Vs)
        # print("Vs_prior: %s" % Vs_prior)
        Bs, forces, min_ns = calculate_bayes_factors(
            zeta_ts=zeta_ts, zeta_sps=zeta_sps, ns=ns, Vs=Vs, Vs_pi=Vs_prior)
        _, forces_grad_only, _ = calculate_bayes_factors(
            zeta_ts=zeta_ts * 0.0, zeta_sps=zeta_sps, ns=ns, Vs=Vs, Vs_pi=Vs_prior)
        _, forces_drift_only, _ = calculate_bayes_factors(
            zeta_ts=zeta_ts, zeta_sps=zeta_sps * 0.0, ns=ns, Vs=Vs, Vs_pi=Vs_prior)

        # # Calculate the Bs that we would get with sufficient n
        # suf_Bs, _, _ = calculate_bayes_factors(zeta_ts = zeta_ts, zeta_sps = zeta_sps, ns = min_ns, Vs = Vs, Vs_pi = Vs_prior)
        # print (suf_Bs)

        # # Check the inferred gradient compared to the simulated gradient

        # Set non-calculated values to nan
        zeta_sps[np.abs(zeta_sps) < abs_tol] = np.nan
        # print(zeta_sps)
        gdt_inf = zeta_sps * vec_to_2D(np.sqrt(Vs))
        median_gdt_inf = np.nanmedian(np.abs(gdt_inf), axis=0)
        gdt_sim = np.asarray([D_0 * k, 0]) * dt

        # Check total force
        alpha_dt_inf = zeta_ts * vec_to_2D(np.sqrt(Vs))
        mean_alpha_dt_inf = np.median(alpha_dt_inf, axis=0)

        # Simulated SNRs
        # alpha_dt_sim = np.asarray([ksi * D_0 * k , 0]) * dt
        # zeta_ts_exp = alpha_dt_sim / vec_to_2D(np.sqrt(Vs))
        zeta_sps_exp = gdt_sim / vec_to_2D(np.sqrt(Vs))

        # Expected order of magnitude of the variance
        # V_magnitude = 4 * D_0 * dt

        # Estimate diffusivity
        D_est = Vs / 4 / dt
        median_D_est = np.median(D_est)

        # General parameters
        # print(np.log10(Bs))
        # print(min_ns)
        # print(min_ns - ns)
        # print((min_ns - ns).astype(int))
        # print(forces)
        print("\n\nMean jump along x: <dx>=\t%.2g um" %
              (np.mean(dr_mean, axis=0)[0]))
        # print ("Expected mean jump along x: <dx>=\t%.2g um" % (ksi * D_0 * k * dt))

        print("\n\nMedian D:\t%.2g um^2/s" % (median_D_est))
        print("Measured range of D:\t[%.2g; %.2g] um^2/s" %
              (np.min(D_est), np.max(D_est)))
        # print ("Expected D range:\t[%.2g; %.2g] um^2/s" % (D_0, D_0 * D_ratio))

        # Inferred gradient
        print("\n\nInferred <|D'|>:\t%s um/s" % (median_gdt_inf / dt))
        # print("Simulated <|D'|>:\t%s um/s" % (gdt_sim / dt))

        print("\n\nMedian zeta_ts:\t%s" % (np.median(zeta_ts, axis=0)))
        # print ("Expected zeta_ts:\t%s" % (np.median(zeta_ts_exp, axis = 0)))
        print("Median abs(zeta_sps):\t%s" %
              (np.nanmedian(np.abs(zeta_sps), axis=0)))
        print("Expected zeta_sps:\t%s" % (np.median(zeta_sps_exp, axis=0)))
        print("Median n:\t%i" % np.median(ns, axis=0))
        print("Median min_n:\t%i" % np.median(min_ns, axis=0))
        print("max(min_n):\t%i" % np.max(min_ns, axis=0))

        # Check total force
        print("\n\nInferred <alpha dt>:\t%s um" % (mean_alpha_dt_inf))
        # print("Simulated <alpha dt>:\t%s um" % (alpha_dt_sim))

        print("\n\n")

        # Prepare to output into a file
        # Output: ksi - first line and next all log10(Bs)
        cells_number = np.size(ns, 0)
        # output = np.zeros((cells_number + 1, 1), dtype = np.float16)
        # output[0] = ksi
        # output[1:] = np.log10(Bs)
        # print(output)
        # print(ns)
        # print(np.log10(Bs))
        output_df = pd.DataFrame(columns=["ksi", "log10_B", "force_evidence", "zeta_t_x", "zeta_t_y",
                                          "zeta_sp_x", "zeta_sp_y", "n_mean", "min_n"], dtype=np.float16)
        output_df["log10_B"] = np.log10(Bs)[:, 0]
        output_df["force_evidence"] = forces[:, 0]
        output_df["n_mean"] = ns[:, 0]
        output_df["min_n"] = min_ns.astype(int)[:, 0]
        output_df["ksi"].loc[0] = ksi
        output_df["zeta_t_x"] = zeta_ts[:, 0]
        output_df["zeta_t_y"] = zeta_ts[:, 1]
        output_df["zeta_sp_x"] = zeta_sps[:, 0]
        output_df["zeta_sp_y"] = zeta_sps[:, 1]
        # output_df = pd.DataFrame(data = {"log10_B": np.log10(Bs)[:, 0],
        # 	"n_mean": ns[:, 0]})
        # output_df.assign(ksi = np.nan)
        # print(output_df)
        # print(output_df.dtypes)

        # Save the file
        dat_file = results_folder + filename + "_" + mesh + '.dat'
        output_df.to_csv(dat_file)
        # with open(dat_file, 'w') as f:
        # 	csv_writer = csv.writer(f, delimiter = CSV_DELIMITER, lineterminator = '\n')
        # 	csv_writer.writerows(output.tolist())

        # Plot
        if bl_output_map:
            cells = analysis_tree[mesh].data  # `cells` contains the mesh

            def png_name(str):
                return results_folder + filename + "_" + mesh + str + '.png'

            # Detected forces
            my_map = pd.DataFrame(forces, index=n.index, columns=[
                                  'Evidence for models'])
            map_plot(my_map, cells=cells,
                     output_file=png_name("_forces"), clip=False)

            my_map = pd.DataFrame(forces_grad_only, index=n.index, columns=[
                                  'Evidence for models | alpha = 0'])
            map_plot(my_map, cells=cells, output_file=png_name(
                "_forces_grad_only"), clip=False)

            my_map = pd.DataFrame(forces_drift_only, index=n.index, columns=[
                                  'Evidence for models | g = 0'])
            map_plot(my_map, cells=cells, output_file=png_name(
                "_forces_drift_only"), clip=False)

            # Log10(B)
            my_map = pd.DataFrame(np.log10(Bs), index=n.index, columns=[
                                  'log10(B) clipped'])
            map_plot(my_map, cells=cells,
                     output_file=png_name("_log_B"), clip=True)

            # D
            my_map = pd.DataFrame(D.T[0], index=n.index, columns=['D'])
            map_plot(my_map, cells=cells,
                     output_file=png_name("_D"), clip=False)

            # Alpha dt
            my_map = pd.DataFrame(alpha_dt_inf, index=n.index, columns=[
                                  '$alpha dt$ x', '$alpha dt$ y'])
            map_plot(my_map, cells=cells, output_file=png_name(
                "_alpha_dt"), clip=False)

            # g dt
            my_map = pd.DataFrame(gdt_inf, index=n.index, columns=[
                                  '$g dt$ x', '$g dt$ y'])
            map_plot(my_map, cells=cells,
                     output_file=png_name("_g_dt"), clip=False)

            # n
            my_map = pd.DataFrame(ns, index=n.index, columns=['n'])
            map_plot(my_map, cells=cells,
                     output_file=png_name("_n"), clip=False)

            #
            # my_map = pd.DataFrame(D.T[0], index = n.index, columns = ['D'])
            # my_map = pd.DataFrame(alpha_dt_inf, index = n.index, columns = ['$alpha dt$ x', '$alpha dt$ y'])
            #

        # # map_plot(my_map, cells=cells, show=False)
        # # ...
        # # plt.show() # waits for the user to close the resulting window
        # ## ... or alternatively plot in a file
        #

        # save the intermediate snr-related maps
        if anything_new:
            save_rwa(rwa_file, analysis_tree, force=True)

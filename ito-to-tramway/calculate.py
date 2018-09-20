#!/usr/bin/env python

"""
Calculate Bayes factors after TRamWAY basic inference if finished.
Plots figures.
Special care is taken of empty cells that disappear from the output of the `infer` procedure.
"""

import csv
import time
from tramway.core import *
from tramway.plot.mesh import plot_voronoi
from tramway.helper import *
import os.path
import numpy as np
import pandas as pd
import matplotlib
# matplotlib.use('Agg')  # enable for console runs with no displays
import matplotlib.pyplot as plt
from set_figure_size import set_figure_size

from bayes_factors.calculate_bayes_factors import calculate_bayes_factors
from constants import *

# for theoretical checks
from constants import D_0, k


def calculate(csv_file, results_folder, bl_produce_maps, dt, snr_label, localization_error, verbose=False):
    # Define the discrete colormap for the Bayes factor
    alpha = 1.0
    col_yes = np.asarray([194,    216, 52]) / 255.0
    col_no = np.asarray([141,  24, 26]) / 255.0
    col_idk = np.asarray([255,   255,    255]) / 255.0
    # colors = np.transpose(np.stack([np.append(red, alpha), np.append(
    #     yellow, alpha), np.append(green, alpha)], axis=1))
    colors = np.transpose(np.stack([col_no, col_idk, col_yes], axis=1))
    # print(colors)
    cmap_bayes_factor = matplotlib.colors.ListedColormap(colors, "bayes_factor")

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

    def sum_cells(a): return np.nansum(a, axis=0, keepdims=True)

    def sum_dims(a): return np.nansum(a, axis=1, keepdims=True)

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
    # analysis_tree = load_rwa(rwa_file)
    analysis_tree = RWAStore(rwa_file, 'r').peek('analyses')
    print(analysis_tree)

    # loop over the available meshes
    anything_new = False
    # `mesh` is a label (a key in dict-like `analysis_tree`)
    for mesh in analysis_tree:  # ['kmeans_20']: - shortcut for the vlp plots
        print("Using ", rwa_file)
        # print("Meshes found: ", analysis_tree[mesh])
        # return

        # infer snr-related maps with the 'snr' plugin
        # (will skip to the next step for precomputed maps)
        if snr_label not in analysis_tree[mesh]:
            # `infer` adds a subtree identified by `snr_label`
            infer(analysis_tree, 'snr', input_label=mesh, output_label=snr_label,
                  max_iter=50, localization_error=localization_error)
            anything_new = True
            save_rwa(rwa_file, analysis_tree, force=True)

        snr = analysis_tree[mesh][snr_label].data

        # get cell centers
        cell_centers = analysis_tree[mesh].data.tessellation.cell_centers  # [snr.maps.index]

        cells_total_len = analysis_tree[mesh].data.location_count.size
        cells_calculated = np.asarray(snr['n'].index)

        # Restore the size of the variables to the original cells number
        n = np.ones([cells_total_len, 1]) * np.nan
        D = np.ones([cells_total_len, 1]) * np.nan
        zeta_total = np.ones([cells_total_len, 2]) * np.nan
        zeta_spurious = np.ones([cells_total_len, 2]) * np.nan
        dr = np.ones([cells_total_len, 2]) * np.nan
        dr2 = np.ones([cells_total_len, 2]) * np.nan

        n[cells_calculated, :] = snr['n'].values
        D[cells_calculated, :] = snr['diffusivity'].values
        zeta_total[cells_calculated, :] = snr['zeta_total'].values
        zeta_spurious[cells_calculated, :] = snr['zeta_spurious'].values
        dr[cells_calculated, :] = snr['dr'].values
        dr2[cells_calculated, :] = snr['dr2'].values

        # print(dr2)
        if True:  # Vpi_name not in snr.variables:
            # `V_prior` is not computed directly by the 'snr' plugin
            # because the plugin's main routine may be independently applied to
            # subsets of cells instead of all the cells at a time
            n_prior = np.nansum(n) - n
            dr_prior = sum_cells(dr) - dr
            dr2_prior = sum_cells(dr2) - dr2

            # calculate biased varaince in current bin
            dr_mean = dr / n
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

        # print(zeta_ts)
        # print(zeta_sps)
        # print(zeta_ts / zeta_sps)
        # Vs =

        # print("Vs: %s" % Vs)
        # print("Vs_prior: %s" % Vs_prior)
        Bs, forces, min_ns = calculate_bayes_factors(
            zeta_ts=zeta_ts, zeta_sps=zeta_sps, ns=ns, Vs=Vs, Vs_pi=Vs_prior)
        # _, forces_grad_only, _ = calculate_bayes_factors(
        #     zeta_ts=zeta_ts * 0.0, zeta_sps=zeta_sps, ns=ns, Vs=Vs, Vs_pi=Vs_prior)
        # _, forces_drift_only, _ = calculate_bayes_factors(
        #     zeta_ts=zeta_ts, zeta_sps=zeta_sps * 0.0, ns=ns, Vs=Vs, Vs_pi=Vs_prior)

        # # Calculate the Bs that we would get with sufficient n
        # suf_Bs, _, _ = calculate_bayes_factors(zeta_ts = zeta_ts, zeta_sps = zeta_sps, ns = min_ns, Vs = Vs, Vs_pi = Vs_prior)
        # print (suf_Bs)

        # # Check the inferred gradient compared to the simulated gradient

        # Set non-calculated values to nan
        zeta_sps[np.abs(zeta_sps) < abs_tol] = np.nan
        # print(zeta_sps)
        gdt_inf = zeta_sps * vec_to_2D(np.sqrt(Vs))
        median_gdt_inf = np.nanmedian(np.abs(gdt_inf), axis=0)
        gdt_abs = np.sqrt(sum_dims(gdt_inf**2))
        gdt_sim = np.asarray([D_0 * k, 0]) * dt

        # Check total force (why inf? why?)
        alpha_dt_inf = zeta_ts * vec_to_2D(np.sqrt(Vs))
        mean_alpha_dt_inf = np.nanmedian(alpha_dt_inf, axis=0)
        # Absolute values
        alpha_dt_abs = np.sqrt(sum_dims(alpha_dt_inf**2))

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
        if verbose:
            print("\n\nMean jump along x: <dx>=\t%.2g um" %
                  (np.nanmean(dr_mean, axis=0)[0]))
            # print ("Expected mean jump along x: <dx>=\t%.2g um" % (ksi * D_0 * k * dt))

            print("\n\nMedian D:\t%.2g um^2/s" % (median_D_est))
            print("Measured range of D:\t[%.2g; %.2g] um^2/s" %
                  (np.nanmin(D_est), np.nanmax(D_est)))
            # print ("Expected D range:\t[%.2g; %.2g] um^2/s" % (D_0, D_0 * D_ratio))

            # Inferred gradient
            print("\n\nInferred <|D'|>:\t%s um/s" % (median_gdt_inf / dt))
            # print("Simulated <|D'|>:\t%s um/s" % (gdt_sim / dt))

            print("\n\nMedian zeta_ts:\t%s" % (np.nanmedian(zeta_ts, axis=0)))
            # print ("Expected zeta_ts:\t%s" % (np.median(zeta_ts_exp, axis = 0)))
            print("Median abs(zeta_sps):\t%s" %
                  (np.nanmedian(np.abs(zeta_sps), axis=0)))
            print("Expected zeta_sps:\t%s" % (np.nanmedian(zeta_sps_exp, axis=0)))
            print("Median n:\t%i" % np.nanmedian(ns, axis=0))
            print("Median min_n:\t%i" % np.nanmedian(min_ns, axis=0))
            print("max(min_n):\t%i" % np.nanmax(min_ns, axis=0))

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
        output_df = pd.DataFrame(columns=["ksi", "x_center", "y_center", "log10_B", "force_evidence", "zeta_t_x",
                                          "zeta_t_y", "zeta_sp_x", "zeta_sp_y", "n_mean", "min_n"], dtype=np.float16)

        output_df["log10_B"] = np.log10(Bs)[:, 0]
        output_df["force_evidence"] = forces[:, 0]
        output_df["n_mean"] = ns[:, 0]
        output_df["min_n"] = min_ns.astype(int)[:, 0]
        output_df["ksi"].loc[0] = ksi
        output_df["zeta_t_x"] = zeta_ts[:, 0]
        output_df["zeta_t_y"] = zeta_ts[:, 1]
        output_df["zeta_sp_x"] = zeta_sps[:, 0]
        output_df["zeta_sp_y"] = zeta_sps[:, 1]
        output_df["x_center"] = cell_centers[:, 0]
        output_df["y_center"] = cell_centers[:, 1]

        # Add the frame also to the .rwa file
        analysis_tree[mesh][snr_label].add(output_df, label='bayes_factors')
        anything_new = True

        # output_df = pd.DataFrame(data = {"log10_B": np.log10(Bs)[:, 0],
        # 	"n_mean": ns[:, 0]})
        # output_df.assign(ksi = np.nan)
        # print(output_df)
        # print(output_df.dtypes)

        # Check that the output folder exists
        if not os.path.isdir(results_folder):
            # Recreate the folder
            try:
                os.makedirs(results_folder)
            except Exception as e:
                print(e)

        # Save the file
        dat_file = os.path.join(results_folder, filename + "_" + mesh + '.dat')
        output_df.to_csv(dat_file)
        # with open(dat_file, 'w') as f:
        # 	csv_writer = csv.writer(f, delimiter = CSV_DELIMITER, lineterminator = '\n')
        # 	csv_writer.writerows(output.tolist())

        # Plot
        if bl_produce_maps:
            cells = analysis_tree[mesh].data  # `cells` contains the mesh

            def png_name(name):

                return os.path.join(results_folder, filename + "_" + mesh + "_" + name + '.png')

            def pdf_name(name):
                return os.path.join(results_folder, filename + "_" + mesh + "_" + name + '.pdf')

            def plot_me(map, name, colormap='inferno', alpha=1.0, bl_plot_mesh=False, mesh_color=[0, 0, 0], colorbar='nice', ticks=False, letter_label=False, colorbar_legend=False, vector=False):
                linewidth = 0.1
                page_width_frac = 1 / 3.0
                pagewidth_in = 6.85
                font_size = 8
                dpi = 100
                alpha_mesh = 0.15

                # the height will be adjusted later
                figsize = np.asarray([3.0, 1.0]) * page_width_frac * pagewidth_in  # in inches
                axis_height = figsize[1]

                figsize = tuple(figsize)

                # Set default figure font size and LaTeX usage
                matplotlib.rcParams.update({'font.size': font_size})

                if not vector:
                    map_plot(map, cells=cells,
                             show=False, clip=False, colormap=colormap, alpha=alpha, linewidth=linewidth, figsize=figsize, dpi=dpi, aspect='equal', colorbar=colorbar)
                else:
                    # linewidth = 0
                    map_plot(map, cells=cells,
                             show=False, clip=False, colormap=colormap, alpha=alpha, markerlinewidth=linewidth, figsize=figsize, dpi=dpi, aspect='equal', colorbar=colorbar, transform=None)

                if bl_plot_mesh:
                    color = tuple(mesh_color + [alpha_mesh])
                    plot_voronoi(cells=cells, color=color,
                                 centroid_style=None, linewidth=linewidth)

                # Manual figure adjustments
                fig = plt.gcf()
                fig.set_dpi(dpi)
                fig.set_figwidth(figsize[0])
                fig.set_figheight(figsize[1])

                # # Enforce a certain axis height
                # ax = fig.axes[0]
                # axis_size = ax.get_position()
                # ax.set_position([axis_size.x0, axis_size.y0, axis_size.width, axis_height])

                # Remove ticks
                if not ticks:
                    plt.tick_params(
                        axis='x',          # changes apply to the x-axis
                        which='both',      # both major and minor ticks are affected
                        bottom=False,      # ticks along the bottom edge are off
                        top=False,         # ticks along the top edge are off
                        labelbottom=False)  # labels along the bottom edge are off
                    plt.tick_params(
                        axis='y',          # changes apply to the x-axis
                        which='both',      # both major and minor ticks are affected
                        left=False,      # ticks along the bottom edge are off
                        right=False,         # ticks along the top edge are off
                        labelleft=False)  # labels along the bottom edge are off

                # add label
                if letter_label:
                    label_location = [0.025, 1.03]
                    # str_label = chr(ord('a') + plot_me.count)
                    # plot_me.count += 1
                    ax = fig.gca()
                    ax.text(label_location[0], label_location[1],
                            letter_label, transform=ax.transAxes, fontsize=font_size)

                # Colorbar legend
                if colorbar_legend:
                    print(fig.axes)
                    ax = fig.axes[1].set_ylabel(colorbar_legend, rotation=90)

                # fig.tight_layout()
                plt.savefig(pdf_name(name), bbox_inches='tight', pad_inches=0)
                plt.savefig(png_name(name), bbox_inches='tight', pad_inches=0)
                fig.clf()

                # map_plot(map, cells=cells,
                #          output_file=pdf_name(name), clip=False, colormap=colormap, alpha=alpha, linewidth=linewidth)
                # colormap = 'inferno'
                # map_plot(map, cells=cells,
                #          output_file=png_name(name), clip=False, colormap='viridis')

            # plot_me.count = 0

            # Alpha dt absolute values
            my_map = pd.DataFrame(alpha_dt_abs, index=n.index, columns=['$\\alpha \Delta t$'])
            plot_me(my_map, "alpha_dt", letter_label='a', colorbar_legend='$\\mu \\mathrm{m}$')

            # Alpha dt vector
            my_map = pd.DataFrame(alpha_dt_inf, index=n.index,  columns=[
                '$\\alpha \Delta t$ x', '$\\alpha \Delta t$ y'])
            plot_me(my_map, "alpha_dt_vec", colormap='inferno', ticks=False,
                    vector=True, colorbar_legend='$\\mu \\mathrm{m}$', bl_plot_mesh=True, mesh_color=[1, 1, 1])

            # Log10(B)
            my_map = pd.DataFrame(np.log10(Bs), index=n.index, columns=[
                '$\log_{10}(B)$'])
            plot_me(my_map, "log_B", letter_label='e')

            # Detected forces
            my_map = pd.DataFrame(forces, index=n.index, columns=['Active force'])
            plot_me(my_map, "bayes_factor", colormap=cmap_bayes_factor,
                    alpha=alpha, bl_plot_mesh=True, colorbar=False, letter_label='f')

            # g dt
            my_map = pd.DataFrame(gdt_abs, index=n.index, columns=['$g \Delta t$'])
            plot_me(my_map, "g_dt", letter_label='d', colorbar_legend='$\\mu \\mathrm{m}$')

            # n
            my_map = pd.DataFrame(ns, index=n.index, columns=['$n$'])
            plot_me(my_map, "n", letter_label='b')

            # D
            my_map = pd.DataFrame(D.T[0], index=n.index, columns=['$D$'])
            plot_me(my_map, "D", letter_label='c', colorbar_legend='$\\mu \\mathrm{m^2/s}$')

            #

            #

            #

            #

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
            print(analysis_tree)

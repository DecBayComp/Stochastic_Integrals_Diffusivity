#!/usr/bin/env python

"""
Calculate Bayes factors after TRamWAY basic inference if finished.
Plot figures.
Special care is taken of empty cells that disappear from the output of the `infer` procedure.
"""

# Special import
if 1:
    import matplotlib
    # matplotlib.use('Agg')  # enable for console runs with no displays

import copy
import csv
import logging
import os.path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm import tqdm

# for theoretical checks
from constants import CSV_DELIMITER, D_0, abs_tol, colors, k
from load_trajectory import load_trajectory
from Sashas import Sashas
from tramway.core.analyses.lazy import Analyses
# from tramway.core import RWAStore
# from tramway.core import save_rwa
from tramway.helper import (RWAStore, cell_plot, infer, infer0, load_rwa,
                            map_plot, save_rwa, tessellate)
# from bayes_factors.calculate_bayes_factors import calculate_bayes_factors
from tramway.inference import Translocations, bayes_factors, distributed
from tramway.plot.mesh import plot_voronoi

cells = []
index = []


def my_save(rwa_file, analysis_tree, **kwargs):
    from tramway.helper import save_rwa
    kwargs['force'] = True
    kwargs['compress'] = False
    save_rwa(rwa_file, analysis_tree, **kwargs)


save_rwa = my_save


def calculate(csv_file, results_folder, bl_produce_maps, snr_label, sigma, tessellation_parameters=False,  verbose=False, recalculate=False, ticks=False, page_width_frac=0.333, min_diffusivity=1e-8):
    # Define the discrete colormap for the Bayes factor
    # alpha = 1.0
    # localization_error is now called sigma
    bayes_factor_label = 'bayes_factor'
    localization_error_BF = sigma**2

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
    except Exception:
        print(".csv file not found. Skipping reading ksi")

    # snr_label = 'snr(mu=0)'
    # Vpi_name = 'V_prior'

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
    analysis_tree = load_trajectory(rwa_file)

    if not analysis_tree:
        # use for vlps
        count = 20
        label = 'kmeans_n_{count}'.format(count=count)
        print(analysis_tree)
        logging.warn(
            'No tessellation found. Performing a standard {label} tessellation'.format(label=label))
        tessellate(analysis_tree, 'kmeans', ref_distance=0, avg_location_count=count,
                   prune=False, label=label)

        save_rwa(rwa_file, analysis_tree, force=True)

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
            infer(analysis_tree, 'd.conj_prior', input_label=mesh, output_label=snr_label,
                  sigma=sigma)
            save_rwa(rwa_file, analysis_tree, force=True)

        snr = analysis_tree[mesh][snr_label].data
        global index
        index = snr['n'].index

        zeta_ts = snr['zeta_total']
        zeta_sps = snr['zeta_spurious']
        ns = snr['n']
        Vs = snr['V']
        D = snr['diffusivity']

        # New Bayes factor calculation procedure
        if recalculate or bayes_factor_label not in analysis_tree[mesh][snr_label]:
            BF_results = infer(analysis_tree, 'bayes_factor',
                               input_label=[mesh, snr_label], output_label=bayes_factor_label, localization_error=localization_error_BF)
            save_rwa(rwa_file, analysis_tree, force=True)
        else:
            BF_results = analysis_tree[mesh][snr_label][bayes_factor_label].data

        lg_Bs, forces, min_ns = [BF_results[label] for label in ['lg_B', 'force', 'min_n']]

        alpha_dt_inf = pd.DataFrame(index=index)
        alpha_dt_inf['x'] = zeta_ts['zeta_total x'] * np.sqrt(Vs.V)
        alpha_dt_inf['y'] = zeta_ts['zeta_total y'] * np.sqrt(Vs.V)
        mean_alpha_dt_inf = np.nanmedian(alpha_dt_inf, axis=0)
        alpha_dt_abs = np.sqrt((alpha_dt_inf**2).sum(axis=1))

        gdt_inf = pd.DataFrame(index=index)
        gdt_inf['x'] = zeta_sps['zeta_spurious x'] * np.sqrt(Vs.V)
        gdt_inf['y'] = zeta_sps['zeta_spurious y'] * np.sqrt(Vs.V)
        gdt_abs = np.sqrt((gdt_inf**2).sum(axis=1))

        # Calculate force masked by the thresholded Bayes factor
        alpha_dt_masked = pd.DataFrame(index=index)
        alpha_dt_masked['x'] = alpha_dt_inf['x'] * (forces.force > 0)
        alpha_dt_masked['y'] = alpha_dt_inf['y'] * (forces.force > 0)

        # Make a simple estimate of lambda where there are no forces
        lambda_est = pd.DataFrame(index=index)
        lambda_est = (zeta_ts.values * zeta_sps.values).sum(axis=1) / (zeta_sps ** 2).sum(axis=1)
        lambda_est[forces.force > 0] = np.nan
        print('Lambda estimates. min: {min:.2f}, median: {median:.2f}, max: {max:.2f}'.format(
            min=np.nanmin(lambda_est),
            median=np.nanmedian(lambda_est),
            max=np.nanmax(lambda_est)))

        # Prepare to output into a file
        # Output: ksi - first line and next all log10(Bs)
        x_centers, y_centers = analysis_tree[mesh].data.tessellation.cell_centers.T

        output_df = pd.DataFrame({"log10_B": lg_Bs['lg_B'],
                                  "force_evidence": forces.force,
                                  "min_n": min_ns['min_n'],
                                  "n_mean": ns.n,
                                  "ksi": ksi,
                                  "zeta_t_x": zeta_ts['zeta_total x'],
                                  "zeta_t_y": zeta_ts['zeta_total y'],
                                  "zeta_sp_x": zeta_sps['zeta_spurious x'],
                                  "zeta_sp_y": zeta_sps['zeta_spurious y'],
                                  'D': D.diffusivity,
                                  'D_CI_lower': snr['ci low']['ci low'],
                                  'D_CI_upper': snr['ci high']['ci high'],
                                  "x_center": x_centers,
                                  "y_center": y_centers,
                                  })

        # Add the data frame also to the .rwa file
        analysis_tree[mesh][snr_label].add(output_df, label='bayes_factors_df')
        anything_new = True

        if results_folder:
            # # Check that the output folder exists
            # if not os.path.isdir(results_folder):
            #     # Recreate the folder
            #     try:
            #         os.makedirs(results_folder)
            #     except Exception as e:
            #         print(e)

            # Save the file
            dat_file = os.path.join(results_folder, filename + "_" + mesh + '.dat')
            output_df.to_csv(dat_file)

            # Plot
            if bl_produce_maps:
                global cells
                cells = analysis_tree[mesh].data  # `cells` contains the mesh

                def fig_name(name):
                    return os.path.join(results_folder, filename + "_" + mesh + "_" + name)

                with tqdm(total=10, desc='Plotting: ') as pbar:

                    # Detected forces
                    plot_me(forces, ['Active force'], fig_name=fig_name(
                        "forces"), letter_label='f', colormap=cmap_bayes_factor, bl_plot_mesh=True, colorbar=False, clim=[-1, 1], ticks=ticks, page_width_frac=page_width_frac)
                    pbar.update()

                    # Alpha dt absolute values
                    plot_me(alpha_dt_abs.to_frame(), ['$\\alpha \Delta t$'],
                            fig_name=fig_name("alpha_dt"), letter_label='c', colormap='inferno',
                            colorbar_legend='$\\mu \\mathrm{m}$', ticks=ticks, page_width_frac=page_width_frac)
                    pbar.update()

                    # Log10(B)
                    plot_me(lg_Bs, ['$\log_{10} K$'], fig_name=fig_name(
                        "lg_B"), letter_label='e', ticks=ticks)
                    pbar.update()

                    # g dt
                    plot_me(gdt_abs.to_frame(), ['$g \Delta t$'], fig_name=fig_name("g_dt"),
                            letter_label=None, colorbar_legend='$\\mu \\mathrm{m}$', ticks=ticks, page_width_frac=page_width_frac)
                    pbar.update()

                    # n
                    plot_me(ns, ['$n$'], fig_name=fig_name("n"), letter_label='b',
                            ticks=ticks, page_width_frac=page_width_frac)
                    pbar.update()

                    # D
                    plot_me(D, ['$D$'], fig_name=fig_name("D"), letter_label='d',
                            colorbar_legend='$\\mu \\mathrm{m^2/s}$', ticks=ticks, page_width_frac=page_width_frac)
                    pbar.update()

                    # # lambda estimate
                    # plot_me(lambda_est.to_frame(), ['$\hat \lambda$'],
                    #         fig_name=fig_name("lambda_est"), ticks=ticks)
                    # pbar.update()

                    # Alpha dt masked
                    plot_me(alpha_dt_masked, [
                        '$\\alpha \Delta t$ x', '$\\alpha \Delta t$ y'],                     fig_name=fig_name("alpha_dt_masked"), colormap='inferno', vector=True, colorbar_legend='$\\mu \\mathrm{m}$', bl_plot_mesh=True, mesh_color=[1, 1, 1], ticks=ticks, page_width_frac=page_width_frac)
                    pbar.update()

                    # Alpha dt vector
                    plot_me(alpha_dt_inf, [
                        '$\\alpha \Delta t$ x', '$\\alpha \Delta t$ y'],                     fig_name=fig_name("alpha_dt_vec"), letter_label='f', colormap='inferno', vector=True, colorbar_legend='$\\mu \\mathrm{m}$', bl_plot_mesh=True, mesh_color=[1, 1, 1], ticks=ticks, page_width_frac=page_width_frac)
                    pbar.update()

                    # g dt vector
                    plot_me(gdt_inf, [
                        '$g \Delta t$ x', '$g \Delta t$ y'],                     fig_name=fig_name("g_dt_vec"), colormap='inferno', vector=True, colorbar_legend='$\\mu \\mathrm{m}$', bl_plot_mesh=True, mesh_color=[1, 1, 1], ticks=ticks, page_width_frac=page_width_frac)
                    pbar.update()

                    # Raw trajectories
                    plot_me(analysis_tree, columns=None, fig_name=fig_name(
                        "raw_trajectories"), letter_label='a', colormap='inferno', bl_plot_mesh=False, colorbar=True, ticks=ticks, max_traj=1000, xlims=[0, 2], ylims=[0, 2], page_width_frac=page_width_frac)
                    pbar.update()

        # save the intermediate snr-related maps
        if anything_new:
            save_rwa(rwa_file, analysis_tree, force=True)
            # print(analysis_tree)


def plot_me(map, columns, fig_name, colormap='inferno', alpha=1.0, bl_plot_mesh=False, mesh_color=[0, 0, 0], colorbar='nice', clim=None, ticks=False, letter_label=False, colorbar_legend=False, vector=False, max_traj=None, xlims=None, ylims=None, page_width_frac=0.333):
    """
    To plot raw trajectories instead of an inferred map, set input to an analysis tree
    """
    linewidth = 0.1

    pagewidth_in = 6.85
    font_size = 8
    dpi = 100
    alpha_mesh = 0.15

    if isinstance(map, np.ndarray):
        map = pd.DataFrame(map, index=index, columns=columns)
    elif isinstance(map, pd.DataFrame) or isinstance(map, pd.Series):
        map = map
        map.columns = columns
    elif isinstance(map, Analyses):
        pass
    else:
        raise TypeError("Unrecoginized map type")

    # the height will be adjusted later
    figsize = np.asarray([3.0, 1.0]) * page_width_frac * pagewidth_in  # in inches
    # axis_height = figsize[1]

    figsize = tuple(figsize)

    # Set default figure font size and LaTeX usage
    matplotlib.rcParams.update({'font.size': font_size})

    if isinstance(map, Analyses):
        plot_raw_trajectories(analysis_tree=map, cmap=colormap,
                              max_traj=max_traj, xlims=xlims, ylims=ylims)
    elif not vector:
        map_plot(map, cells=cells,
                 show=False, clip=False, colormap=colormap, alpha=alpha, linewidth=linewidth, figsize=figsize, dpi=dpi, aspect='equal', colorbar=colorbar, clim=clim)
    else:
        # linewidth = 0
        map_plot(map, cells=cells,
                 show=False, clip=False, colormap=colormap, alpha=alpha, markerlinewidth=linewidth, figsize=figsize, dpi=dpi, aspect='equal', colorbar=colorbar, clim=clim, transform=None)

    if bl_plot_mesh:
        color = tuple(mesh_color + [alpha_mesh])
        plot_voronoi(cells=cells, color=color,
                     centroid_style=None, linewidth=linewidth)

    # auto = plt.gcf()
    # print(auto.number, id(auto))
    # auto.savefig('test4.png')
    # print(plt.get_fignums())
    # print(plt.figure(1).number, id(plt.figure(1)))
    # plt.figure(1).savefig('test5.png')

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
            axis='x',           # changes apply to the x-axis
            which='both',       # both major and minor ticks are affected
            bottom=False,       # ticks along the bottom edge are off
            top=False,          # ticks along the top edge are off
            labelbottom=False)  # labels along the bottom edge are off
        plt.tick_params(
            axis='y',           # changes apply to the y-axis
            which='both',       # both major and minor ticks are affected
            left=False,         # ticks along the bottom edge are off
            right=False,        # ticks along the top edge are off
            labelleft=False)    # labels along the bottom edge are off
    else:
        plt.xlabel('$x, \mu\mathrm{m}$')
        plt.ylabel('$y, \mu\mathrm{m}$')

    # add label
    if letter_label:
        label_location = [0.025, 1.03]
        # str_label = chr(ord('a') + plot_me.count)
        # plot_me.count += 1
        ax = fig.gca()
        ax.text(label_location[0], label_location[1],
                letter_label, transform=ax.transAxes, fontsize=font_size)

    # Colorbar legend
    if colorbar and colorbar_legend:
        # print(fig.axes)
        # ax_bkp = ax
        fig.axes[1].set_ylabel(colorbar_legend, rotation=90)
    # plt.gcf().savefig('test6.png')
    # print(id(plt.gcf()))
    # print(id(fig))

    # fig.tight_layout()

    # # For JBs illustration, add a lign for D
    # ax = fig.gca()
    # ax.plot([0, 1], [0.5] * 2, 'w-', lw=0.5)

    try:
        fig.savefig(fig_name + '.pdf', bbox_inches='tight', pad_inches=0)
    except:
        logging.warn('Unable to save figure. The file might be open')

    # Save zoomed .png version
    factor = 5
    fig.set_figwidth(figsize[0] * factor)
    fig.set_figheight(figsize[1] * factor)
    try:
        fig.savefig(fig_name + '.png', pad_inches=0, bbox_inches='tight')  # )
    except:
        logging.warn('Unable to save figure. The file might be open')
    # plt.close(fig)

    # map_plot(map, cells=cells,
    #          output_file=pdf_name(name), clip=False, colormap=colormap, alpha=alpha, linewidth=linewidth)
    # colormap = 'inferno'
    # map_plot(map, cells=cells,
    #          output_file=png_name(name), clip=False, colormap='viridis')


def plot_raw_trajectories(analysis_tree, cmap, max_traj=None, xlims=None, ylims=None):
    data = analysis_tree.data
    Ni = int(data.n.max())
    # print(data)
    # max_traj = 1200
    line_width = 0.5
    alpha = 1
    cmap = 'inferno'
    xlim = [0, 2]
    ylim = [0, 2]

    # Time and color mesh
    t_min = data.t.min()
    t_max = data.t.max()
    dt = data.iloc[0]['dt']
    t = np.arange(t_min, t_max + dt, dt)
    t_segment_centers = (t[1:] + t[:-1]) / 2
    Nt = len(t)
    color = (t[1:] + t[:-1]) / 2 / t_max

    # Subsample to unclutter
    np.random.seed(0)
    ns = list(range(0, Ni + 1))
    # print(ns)
    np.random.shuffle(ns)
    if max_traj is not None and Ni - 1 > max_traj:
        ns = ns[0:max_traj]
    Ni = len(ns)

    # print(Ni)
    # Create an array of all points of each trajectory including nan
    points = np.zeros((Ni, Nt, 2)) * np.nan
    for i, n in enumerate(ns):
        slice = copy.deepcopy(data[data.n == n].loc[:, ['t', 'x', 'y', 'dt', 'dx', 'dy']])
        # Add an endpoint
        slice = slice.append(slice.iloc[-1][['t', 'x', 'y']] +
                             slice.iloc[-1][['dt', 'dx', 'dy']].values)
        inds = np.rint(slice.t.values / dt).astype(np.int)
        points[i, inds, 0] = slice.x
        points[i, inds, 1] = slice.y
    # print(points)

    def make_segments(points):
        x = points[:, 0]
        y = points[:, 1]
        # Add axis 1
        segments = np.vstack([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([segments[:-1], segments[1:]], axis=1)
        # Return size is N x 2 x 2, where axis == 1 corresponds to the length of the segment. Axis = 2 is the x, y components
        return segments

    # Convert points to line segments
    segments_all = [make_segments(item) for item in points]
    segments_all = np.vstack(segments_all)
    colors_all = np.tile(color, Ni)
    # Drop Nan segments
    nan_inds = np.any(np.isnan(segments_all), (2, 1))
    segments = segments_all[~nan_inds, :, :]
    colors = colors_all[~nan_inds]

    # Plot
    fig, ax = plt.subplots(num=1, clear=True)
    ax.set_aspect('equal')
    coll = LineCollection(segments, cmap=cmap, lw=line_width, alpha=alpha)
    coll.set_array(colors)
    ax.add_collection(coll)

    # ax.set_facecolor('k')

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=t_min, vmax=t_max))
    sm._A = []
    # plt.colorbar(sm)
    colorbar_legend = 't, s'
    try:
        gca_bkp = plt.gca()
        divider = make_axes_locatable(fig.gca())
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(sm, cax=cax)
        cax.set_ylabel(colorbar_legend, rotation=90)
        plt.sca(gca_bkp)
    except AttributeError as e:
        logging.warning(e.args[0], RuntimeWarning)

    if xlims is None:
        xlims = (np.nanmin(points[:, :, 0]), np.nanmax(points[:, :, 0]))
    if ylims is None:
        ylims = (np.nanmin(points[:, :, 1]), np.nanmax(points[:, :, 1]))
    plt.ylim(ylims)
    plt.xlim(xlims)

    plt.title('{} trajectories'.format(Ni))

    return None

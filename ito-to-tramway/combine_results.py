"""
Combine dat files with individual Bayes factor calculations and average the statistics over the trials
"""

# from constants import k, D_0, D_ratio, dt
import glob
import os

import numpy as np
import pandas as pd
from tqdm import tqdm, trange  # console progress meter

from constants import D_0, D_ratio, L
from get_expected_B import get_expected_B

# from tramway.inference.bayes_factors import calculate_bayes_factors


def combine_results(bl_force_reload=False):
    # Constants
    # bl_force_reload=True
    lg_B_abs_treshold = 1.0
    ksi_precision = 2
    quantile = 0.05
    combined_data_filename = "combined_data.dat"
    lgB_filename = "combined_lgB.dat"
    folder_no_perp = r'd:\calculated_data\sim_performance_2D_no_perp'
    folder_with_perp = r'd:\calculated_data\sim_performance_2D_with_perp'
    # folder_with_perp = r'd:\calculated_data\sim_performance_2D_with_perp_1000_internal_steps'
    folders = [folder_no_perp, folder_with_perp]

    # initialize
    avg_data = []
    top_quantile = []
    bottom_quantile = []
    ksis_unique = []
    trials_number = []
    expect_mean_n = []

    # data_with_perp = pd.DataFrame()
    # data_no_perp = pd.DataFrame()
    data = []
    data_lgB = []
    for folder_ind in range(len(folders)):
        # folder_ind = 0
        folder = folders[folder_ind]
        # folder = folder_no_perp

        # Get a list of results files
        results_files = glob.glob(os.path.join(folder, "dat/*.dat"))
        stats_file = os.path.join(folder, combined_data_filename)
        lgB_file = os.path.join(folder, lgB_filename)
        print("\nProcessing folder: ", folder)
        files_count = len(results_files)
        files_count = np.min([files_count, 101 * 100])  # take max trials 101 * 1000

        max_bin = estimate_max_bin_number(results_files)

        # Initialize
        left_ksis = np.zeros(files_count, dtype=np.float32) * np.nan
        true_ksis = np.zeros((files_count, 2), dtype=np.float32) * np.nan
        mean_ns = np.zeros((files_count, 2), dtype=np.float32) * np.nan
        zsp_x_mean = np.zeros((files_count, 2), dtype=np.float32) * np.nan
        zt_y_mean = np.zeros((files_count, 2), dtype=np.float32) * np.nan

        # Dimensions of tr_B: file, cell half, B region
        tr_Bs = np.zeros((files_count, 2, 3), dtype=np.float32) * np.nan
        lg_Bs = pd.DataFrame(columns=['ksi', 'lg_B'], index=np.arange(
            max_bin * files_count), dtype=np.float)

        # Load data
        print('Found %i files' % files_count)
        if bl_force_reload or not os.path.exists(stats_file) or not os.path.exists(lgB_file):
            j = 0
            for i in trange(files_count):
                # i=0
                file = results_files[i]
                df = pd.read_csv(file)
                bins = len(df)

                # try:
                # calculate separately for the left and right half
                x_centers = df["x_center"]
                left_ksis[i] = df["ksi"].loc[0]
                true_ksis[i, :] = np.asarray([1, -1]) * left_ksis[i]
                # halves indicator has True in column 0 if the bin is in the left half
                halves_indicators = np.transpose(np.stack((np.asarray(
                    x_centers) <= 0.5, np.asarray(x_centers) > 0.5)))
                for half in range(2):
                    half_cells = df.loc[halves_indicators[:, half]]

                    # Calculate
                    # ksis in the left half of the box. In the right half they have an opposite sign
                    mean_ns[i, half] = np.nanmean(half_cells["n_mean"])
                    zt_y_mean[i, half] = np.nanmean(half_cells["zeta_t_y"])
                    zsp_x_mean[i, half] = np.nanmean(half_cells["zeta_sp_x"])

                    # Get the fraction of active force detection
                    if 'log_10_B' in half_cells.columns:
                        lgB_name = 'log_10_B'
                    else:
                        lgB_name = 'log10_B'
                    cur_lg_Bs = np.asarray(half_cells[lgB_name])
                    cells_count = np.sum(~np.isnan(half_cells[lgB_name]))
                    # print(cells_count)

                    # calculate half fractions
                    tr_Bs[i, half, 2] = np.sum(
                        (cur_lg_Bs >= lg_B_abs_treshold) * 1.0) / cells_count
                    tr_Bs[i, half, 0] = np.sum(
                        (cur_lg_Bs <= -lg_B_abs_treshold) * 1.0) / cells_count
                    tr_Bs[i, half, 1] = 1 - tr_Bs[i, half, 0] - tr_Bs[i, half, 2]

                # Save Bayes factors
                index = np.arange(j, j + bins)
                # Bs.loc[index]
                lg_Bs.loc[index, 'lg_B'] = df[lgB_name].values
                lg_Bs.loc[index, 'D_sim'] = D_sim(x_centers.values)
                lg_Bs.loc[index, 'ksi'] = true_ksis[i, 1]
                lg_Bs.loc[index[halves_indicators[:, 0]], 'ksi'] = true_ksis[i, 0]
                lg_Bs.loc[index, 'n'] = df.n_mean.values
                lg_Bs.loc[index, 'bl_zt_perp'] = folder_ind

                # lg_Bs.loc[index, 'lg_B_expected'] = get_expected_B(
                #     lg_Bs.loc[index, 'ksi'].values,
                #     lg_Bs.loc[index, 'D_sim'].values, df.n_mean.values, folder_ind)
                # print(lg_Bs.loc[index, 'lg_B_expected'])
                # return
                j += bins

                # except Exception as e:
                #     print("Warning: encountered error while processing file %s. Skipping.\n\n" % (file))
                #     print(e)

            # Convert to an appropriate data frame
            # ksi_rounded is a rounded value of ksi in these simulations
            # left and right halves will appear as independent entries on this table, with their own ksis

            # Cut the Bs data frame
            lg_Bs = lg_Bs.drop(lg_Bs.index[j - 1:])
            # Bs.drop(Bs.index[j - 1:]).to_csv('temp.csv')

            # add each half
            data.append(pd.DataFrame())
            for half in range(2):
                add_data = pd.DataFrame({
                    "ksi": true_ksis[:, half],
                    "ksi_rounded": np.asarray(true_ksis[:, half] * 10**ksi_precision, dtype=np.int),
                    "mean_n": mean_ns[:, half],
                    "zeta_sp_x_mean": zsp_x_mean[:, half],
                    "zeta_sp_x_abs_mean": np.abs(zsp_x_mean[:, half]),
                    "zeta_t_y_mean": zt_y_mean[:, half],
                    "frac_cons": tr_Bs[:, half, 2],
                    "frac_uncert": tr_Bs[:, half, 1],
                    "frac_spur": tr_Bs[:, half, 0]},
                )

                # count the number of trials from 1 half (otherwise they double)
                if half == 0:
                    trials = add_data.groupby(by="ksi_rounded").count().mean().iloc[0]
                    # print(trials)
                    # add_data['trials'] = np.nan
                    add_data.loc[0, 'trials'] = trials
                    trials_number.append(np.mean(trials))

                data[folder_ind] = data[folder_ind].append(add_data, sort=False)

            # reset row counter
            data[folder_ind] = data[folder_ind].reset_index(drop=True)

            # Calculate the expected Bayes factors
            lg_Bs['lg_B_expected'] = get_expected_B(ksis=lg_Bs.ksi.values,
                                                    D_sims=lg_Bs.D_sim.values,
                                                    ns=lg_Bs.n.values, bl_ztpers=lg_Bs.bl_zt_perp.values)

            # save data to csvs
            data[folder_ind].to_csv(stats_file)
            data_lgB.append(lg_Bs)
            lg_Bs.to_csv(lgB_file)

            # Clean up
            left_ksis = []
            cur_lg_Bs = []
            cells_count = []
            tr_Bs = []
        else:
            # Read lg_B
            data_lgB.append(pd.read_csv(lgB_file))

            # Read other data
            cur_data = pd.read_csv(stats_file)
            data.append(cur_data)
            trials_number.append(cur_data.loc[0, 'trials'])
            cur_data = []
            print("Warning: trajectories not reloaded")

        # Average over trials now treating the halves independently
        # print(data.groupby(by="ksi_rounded").mean())
        avg_data.append(data[folder_ind].groupby(by="ksi_rounded").mean())
        top_quantile.append(data[folder_ind].groupby(by="ksi_rounded").quantile(q=0.975))
        bottom_quantile.append(data[folder_ind].groupby(by="ksi_rounded").quantile(q=0.225))
        ksis_unique.append(avg_data[folder_ind]["ksi"])
        expect_mean_n.append(np.mean(data[folder_ind]["mean_n"]))

        print("Folder processed!")
        print("Files analyzed: %i" % files_count)
        print("Trials processed: %.2f" % trials_number[folder_ind])
    return data, data_lgB, ksis_unique, avg_data, expect_mean_n, trials_number


def get_max_bin_number(results_files):
    max_bins = 0
    for file in tqdm(results_files):
        # file = r'd:\calculated_data\sim_performance_2D_no_perp\dat\sim_data_000000001_gwr.dat'
        with open(file) as f:
            for bins, _ in enumerate(f):
                pass
        max_bins = np.max([bins, max_bins])
    return max_bins


def estimate_max_bin_number(results_files):
    max_bins = 0
    look_at = 100
    factor = 2

    if len(results_files) > look_at:
        results_files = results_files[0:look_at]
    for file in results_files:
        # file = r'd:\calculated_data\sim_performance_2D_no_perp\dat\sim_data_000000001_gwr.dat'
        with open(file) as f:
            for bins, _ in enumerate(f):
                pass
        max_bins = np.max([bins, max_bins])
    return factor * max_bins


def D_sim(x):
    return D_0 * (D_ratio - 2 * (D_ratio - 1) * np.abs(x - L / 2))

"""
Combine dat files with individual Bayes factor calculations and average the statistics over the trials
"""

# from constants import k, D_0, D_ratio, dt
import glob
import numpy as np
import os
import pandas as pd
from tqdm import tqdm, trange   # console progress meter


def combine_results(bl_force_reload=False):
    # Constants
    lg_B_abs_treshold = 1.0
    ksi_precision = 2
    quantile = 0.05
    combined_data_filename = "combined_data.dat"
    folder_no_perp = r'd:\calculated_data\sim_performance_2D_no_perp'
    folder_with_perp = r'd:\calculated_data\sim_performance_2D_with_perp'
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
    for f_ind in range(len(folders)):
        folder = folders[f_ind]

        # Get a list of results files
        results_files = glob.glob(os.path.join(folder, "dat/*.dat"))
        stats_file = os.path.join(folder, combined_data_filename)
        print("\nProcessing folder: ", folder)
        files_count = len(results_files)
        files_count = np.min([files_count, 101 * 10])  # take max trials

        # Initialize
        left_ksis = np.zeros(files_count, dtype=np.float32) * np.nan
        true_ksis = np.zeros((files_count, 2), dtype=np.float32) * np.nan
        mean_ns = np.zeros((files_count, 2), dtype=np.float32) * np.nan
        zsp_x_mean = np.zeros((files_count, 2), dtype=np.float32) * np.nan
        zt_y_mean = np.zeros((files_count, 2), dtype=np.float32) * np.nan

        # Dimensions of tr_B: file, cell half, B region
        tr_Bs = np.zeros((files_count, 2, 3), dtype=np.float32) * np.nan

        # Load data
        if bl_force_reload or not os.path.exists(stats_file):

            for i in trange(files_count):
                file = results_files[i]
                df = pd.read_csv(file)

                try:
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
                        cur_lg_Bs = np.asarray(half_cells["log10_B"])
                        cells_count = np.sum(~np.isnan(half_cells["log10_B"]))
                        # print(cells_count)

                        # calculate half fractions
                        tr_Bs[i, half, 2] = np.sum(
                            (cur_lg_Bs >= lg_B_abs_treshold) * 1.0) / cells_count
                        tr_Bs[i, half, 0] = np.sum(
                            (cur_lg_Bs <= -lg_B_abs_treshold) * 1.0) / cells_count
                        tr_Bs[i, half, 1] = 1 - tr_Bs[i, half, 0] - tr_Bs[i, half, 2]
                except Exception as e:
                    print("Warning: encountered error while processing file %s. Skipping.\n\n" % (file))
                    print(e)

            # Convert to an appropriate data frame
            # ksi_rounded is a rounded value of ksi in these simulations
            # left and right halves will appear as independent entries on this table, with their own ksis

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

                data[f_ind] = data[f_ind].append(add_data)

            # reset row counter
            data[f_ind] = data[f_ind].reset_index(drop=True)

            # save data to csvs
            data[f_ind].to_csv(stats_file)

            # Clean up
            left_ksis = []
            cur_lg_Bs = []
            cells_count = []
            tr_Bs = []
        else:
            cur_data = pd.read_csv(stats_file)
            data.append(cur_data)
            trials_number.append(cur_data.loc[0, 'trials'])
            cur_data = []
            print("Warning: trajectories not reloaded")

        # Average over trials now treating the halves independently
        # print(data.groupby(by="ksi_rounded").mean())
        avg_data.append(data[f_ind].groupby(by="ksi_rounded").mean())
        top_quantile.append(data[f_ind].groupby(by="ksi_rounded").quantile(q=0.975))
        bottom_quantile.append(data[f_ind].groupby(by="ksi_rounded").quantile(q=0.225))
        ksis_unique.append(avg_data[f_ind]["ksi"])
        expect_mean_n.append(np.mean(data[f_ind]["mean_n"]))

        print("Folder processed!")
        print("Files analyzed: %i" % files_count)
        print("Trials processed: %.2f" % trials_number[f_ind])
    return data, ksis_unique, avg_data, expect_mean_n, trials_number

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
        # for folder in folders:
        folder = folders[f_ind]

        # Get a list of results files
        results_files = glob.glob(os.path.join(folder, "dat/*.dat"))
        stats_file = os.path.join(folder, combined_data_filename)
        print("\nProcessing folder: ", folder)
        files_count = len(results_files)
        # files_count = 202

        # Load data
        if bl_force_reload or not os.path.exists(stats_file):
            ksis = np.zeros(files_count, dtype=np.float32)
            mean_ns = np.zeros(files_count, dtype=np.float32)
            zsp_abs_mean = np.zeros(files_count, dtype=np.float32)
            zt_y_mean = np.zeros(files_count, dtype=np.float32)
            tr_Bs = np.zeros((files_count, 3), dtype=np.float32)
            # i=10
            progress = 0
            # progress_bar = FloatProgress(min=0, max=100, description="Calculating:")
            # display(progress_bar)
            for i in trange(files_count):
                file = results_files[i]
                df = pd.read_csv(file)

                # Calculate
                ksis[i] = df["ksi"].loc[0]
                mean_ns[i] = np.nanmean(df["n_mean"])
                zsp_abs_mean[i] = np.nanmean(
                    np.sqrt(df["zeta_sp_x"] ** 2.0 + df["zeta_sp_y"] ** 2.0))
                zt_y_mean[i] = np.nanmean(df["zeta_t_y"])

                # Get the fraction of active force detection
                cur_lg_Bs = np.asarray(df["log10_B"])
                cells_count = np.size(cur_lg_Bs, 0)
                tr_Bs[i, 2] = np.sum(
                    (cur_lg_Bs >= lg_B_abs_treshold) * 1.0) / cells_count
                tr_Bs[i, 0] = np.sum(
                    (cur_lg_Bs <= -lg_B_abs_treshold) * 1.0) / cells_count
                tr_Bs[i, 1] = 1 - tr_Bs[i, 0] - tr_Bs[i, 2]

                # if i / files_count * 100.0 >= progress:
                #     progress += 1
                #     print("Progress: %i%% (%i/%i)" % (progress, i, files_count))
                # progress_bar.value = progress

            # Convert to an appropriate data frame
            # ksi_rounded is a rounded value of ksi in these simulations
            data.append(pd.DataFrame({
                "ksi": ksis,
                "ksi_rounded": np.asarray(ksis * 10**ksi_precision, dtype=np.int),
                "mean_n": mean_ns,
                "zeta_sp_abs_mean": zsp_abs_mean,
                "zeta_t_y_mean": zt_y_mean,
                "frac_cons": tr_Bs[:, 2],
                "frac_uncert": tr_Bs[:, 1],
                "frac_spur": tr_Bs[:, 0]},
            ))

            # save data to csvs
            data[f_ind].to_csv(stats_file)

            # Clean up
            ksis = []
            cur_lg_Bs = []
            cells_count = []
            tr_Bs = []
        else:
            data.append(pd.read_csv(stats_file))
            print("Warning: trajectories not reloaded")

        # Average over trials
        # print(data.groupby(by="ksi_rounded").mean())
        avg_data.append(data[f_ind].groupby(by="ksi_rounded").mean())
        top_quantile.append(data[f_ind].groupby(by="ksi_rounded").quantile(q=0.975))
        bottom_quantile.append(data[f_ind].groupby(by="ksi_rounded").quantile(q=0.225))
        ksis_unique.append(avg_data[f_ind]["ksi"])
        trials_number.append(data[f_ind].groupby(by="ksi_rounded").count().iloc[0].tolist()[0])
        expect_mean_n.append(np.mean(data[f_ind]["mean_n"]))

        print("Folder processed!")
        print("Files analyzed: %i" % files_count)
        print("Trials processed: %i" % trials_number[f_ind])
    return data, ksis_unique, avg_data, expect_mean_n, trials_number

"""
This interactve python file analyzes and combines data for the article statistical performance plot and then makes the plot
"""

# %% Magic functions
%matplotlib
%load_ext autoreload
%autoreload 2


# %% Combine results
from constants import output_folder, k, D_0, D_ratio, dt
from bayes_factors.find_marginalized_zeta_t_roots import find_marginalized_zeta_t_roots
import glob
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from set_figure_size import set_figure_size

# For progress visualization
from ipywidgets import FloatProgress
from IPython.display import display

# def combine_results():

# Constants
lg_B_abs_treshold = 1.0
alpha = 0.13
fill_color_sequence = ["r", "b", "b"]
ksi_precision = 2
quantile = 0.05
combined_data_file = "combined_data.csv"

# Get a list of results files
results_files = glob.glob(output_folder + "*.dat")
files_count = len(results_files)
# files_count = 200

# Load data
if not os.path.exists(combined_data_file):
    ksis = np.zeros(files_count, dtype=np.float16)
    mean_ns = np.zeros(files_count, dtype=np.float16)
    tr_Bs = np.zeros((files_count, 3), dtype=np.float16)
    # i=10
    progress = 0
    progress_bar = FloatProgress(min=0, max=100, description="Calculating:")
    display(progress_bar)
    for i in range(files_count):
        file = results_files[i]
        df = pd.read_csv(file)

        # Get ksi
        ksis[i] = df["ksi"].loc[0]
        mean_ns[i] = np.nanmean(df["n_mean"])

        # Get tr_Bs
        cur_lg_Bs = np.asarray(df["log10_B"])
        cells_count = np.size(cur_lg_Bs, 0)
        tr_Bs[i, 2] = np.sum(
            (cur_lg_Bs >= lg_B_abs_treshold) * 1.0) / cells_count
        tr_Bs[i, 0] = np.sum(
            (cur_lg_Bs <= -lg_B_abs_treshold) * 1.0) / cells_count
        tr_Bs[i, 1] = 1 - tr_Bs[i, 0] - tr_Bs[i, 2]

        if i / files_count * 100.0 >= progress:
            progress += 1
#             print("Progress: %i%% (%i/%i)" % (progress, i, files_count))
        progress_bar.value = progress

    # Convert to a data frame
    data = pd.DataFrame({
        "ksi": ksis,
        "ksi_rnd": np.asarray(ksis * 10**ksi_precision, dtype=np.int),
        "mean_n": mean_ns,
        "frac_cons": tr_Bs[:, 2],
        "frac_uncert": tr_Bs[:, 1],
        "frac_spur": tr_Bs[:, 0]},
    )

    # Clean up
    ksis = []
    cur_lg_Bs = []
    cells_count = []
    tr_Bs = []
else:
    data = pd.read_csv(combined_data_file)
    print("Warning: trajectories not reloaded")

# Average over trials
avg_data = data.groupby(by="ksi_rnd").mean()
top_quantile = data.groupby(by="ksi_rnd").quantile(q=0.975)
bottom_quantile = data.groupby(by="ksi_rnd").quantile(q=0.225)
ksis_unique = avg_data["ksi"]
trials_number = data.groupby(by="ksi_rnd").count().iloc[0].tolist()[0]

print("Finished!")
print("Files analyzed: %i" % files_count)
print("Trials processed: %i" % trials_number)

# Save data to csvs
data.to_csv(combined_data_file)


# print(ksis_unique)

# print(avg_data)


# %% >>> Theoretical estimates <<<
n_pi = 4
u = 1.0
dim = 2
theor_Bs = [0.1, 10.0]
D_max = D_0 * D_ratio
sim_abs_zeta_sp = k * D_0 * np.sqrt(dt) / (np.sqrt(D_0) + np.sqrt(D_max))
exp_mean_n = np.mean(data["mean_n"])


# Zeta_t roots calculation (probably needs some testing)
exp_zeta_ts = np.zeros(4, dtype=np.float16)
exp_zeta_ts_over_zeta_sps = np.zeros(4, dtype=np.float16)
# Low B threshold
zeta_t_perp = 0.0
exp_zeta_ts[[1, 2]] = find_marginalized_zeta_t_roots(zeta_sp=[sim_abs_zeta_sp],
                                                     n=exp_mean_n, n_pi=n_pi, B=theor_Bs[0], u=u, dim=dim, zeta_t_perp=zeta_t_perp)
exp_zeta_ts[[0, 3]] = find_marginalized_zeta_t_roots(zeta_sp=[sim_abs_zeta_sp],
                                                     n=exp_mean_n, n_pi=n_pi, B=theor_Bs[1], u=u, dim=dim, zeta_t_perp=zeta_t_perp)

exp_zeta_ts_over_zeta_sps = exp_zeta_ts / sim_abs_zeta_sp
print(sim_abs_zeta_sp)
print(exp_zeta_ts)


# %% >> > Plot << <


# Try to plot from terminal
fig = set_figure_size(num=1, rows=1, page_width_frac=0.5)
_, ax = plt.subplots(num=1)
ax.plot(ksis_unique.tolist(),
        avg_data["frac_spur"].tolist(), color='r', label="spurious")
ax.plot(ksis_unique.tolist(),
        avg_data["frac_uncert"].tolist(), color='g', label="ins. ev.")
ax.plot(ksis_unique.tolist(),
        avg_data["frac_cons"].tolist(), color='b', label="conserv.")

# Prepare rectangles to fill in the format x, y, width, height
xlims = ax.get_xlim()
ylims = ax.get_ylim()
fill_regions = np.zeros((3, 4), dtype=np.float16)
fill_regions[0, :] = [exp_zeta_ts_over_zeta_sps[1], ylims[0],
                      exp_zeta_ts_over_zeta_sps[2] - exp_zeta_ts_over_zeta_sps[1], ylims[1] - ylims[0]]
fill_regions[1, :] = [xlims[0], ylims[0],
                      exp_zeta_ts_over_zeta_sps[0] - xlims[0], ylims[1] - ylims[0]]
fill_regions[2, :] = [exp_zeta_ts_over_zeta_sps[3], ylims[0],
                      xlims[1] - exp_zeta_ts_over_zeta_sps[3], ylims[1] - ylims[0]]

# ax.legend()

# Add fill
for i in range(3):
    ax.add_patch(matplotlib.patches.Rectangle((fill_regions[i, 0], fill_regions[i, 1]),
                                              fill_regions[i, 2], fill_regions[i,
                                                                               3], facecolor=fill_color_sequence[i],
                                              edgecolor="none", fill=True, alpha=alpha))


# Adjust
ax.set_xlabel("Simulated $\zeta_t / \zeta_{sp}$")
ax.set_ylabel("Fraction")
ax.set_title("$\\bar{n} = %i$, trials = %i" %
             (round(exp_mean_n), trials_number))
ax.legend(loc="lower left")
plt.ion()

fig.tight_layout()
# fig.subplots_adjust(bottom=0.1)
plt.show()


# fig.savefig("simulation_results.png")
# fig.savefig("simulation_results.pdf")

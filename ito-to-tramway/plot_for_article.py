"""
Plot statistical data together with theoretical predictions.
Since the horizontal axis features the abs. value of zeta_sp, the regions are duplicated for the opposite signs.
Maybe I should simplify the presentation, recalculate and show the signed values? This would reduce the number of peaks and simplify presentation.
"""

from constants import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from set_figure_size import set_figure_size


def plot_for_article(ksis_unique, avg_data, exp_zeta_ts_over_zeta_sps, expext_mean_n, trials_number):

    # Constants
    alpha = 0.25
    label_location = [0.025, 0.9]
    font_size = 8
    fill_color_sequence = [green, yellow, red, yellow, green]
    rows = 2
    cols = 1
    ylims = [-0.02, 1.02]
    xlims = np.asarray([-1, 1]) * 25.5
    # xlims = [[-25, 25], [-5, 5]]

    fig = set_figure_size(num=1, rows=rows, page_width_frac=0.5)
    _, ax_arr = plt.subplots(rows, cols, num=1, sharex=False)

    count = 0
    for i in range(rows):
        # print(ax_arr[0])
        ax = ax_arr[i]
        ax.plot(ksis_unique[i].tolist(),
                avg_data[i]["frac_spur"].tolist(), color=red, label="sp")
        ax.plot(ksis_unique[i].tolist(),
                avg_data[i]["frac_uncert"].tolist(), color=yellow, label="insf")
        ax.plot(ksis_unique[i].tolist(),
                avg_data[i]["frac_cons"].tolist(), color=green, label="act")

        # Prepare rectangles to fill in the format x, y, width, height
        ax.set_xlim(xlims)
        # print(xlims)
        # ylims = ax.get_ylim()

        # print(np.shape(exp_zeta_ts_over_zeta_sps))
        # print(np.shape([[0], [0]]))
        regions_x_coords = np.concatenate(
            (np.asarray([[1], [1]]) * xlims[0], exp_zeta_ts_over_zeta_sps, np.asarray([[1], [1]]) * xlims[1]), axis=1)

        # Add fill
        for j in range(5):
            ax.add_patch(matplotlib.patches.Rectangle(xy=(regions_x_coords[i, j], ylims[0]), width=regions_x_coords[i, j + 1] - regions_x_coords[i, j],
                                                      height=ylims[1] - ylims[0], facecolor=fill_color_sequence[j], edgecolor="none", fill=True, alpha=alpha))

        # Adjust
        # Mean over ksi values
        # print(len(avg_data))
        # print(avg_data[i])
        ax.set_ylim(ylims)
        avg_zt_y = np.mean(avg_data[i]["zeta_t_y_mean"])
        avg_abs_zsp_x = np.mean(avg_data[i]["zeta_sp_x_abs_mean"])
        print("Now zt_per/zsp is %.3f" % (avg_zt_y / avg_abs_zsp_x))
        print("Advice: set your zt_par/|zsp| ratio to %.3f if you wish to get zt_per = 0.1" %
              (0.1 / avg_abs_zsp_x))

        # Add a label to each plot
        str_label = chr(ord('a') + count)
        ax.text(label_location[0], label_location[1],
                str_label, transform=ax.transAxes, fontsize=font_size)
        count += 1

        # round zt_y
        if np.abs(avg_zt_y) < 1e-3:
            avg_zt_y = 0

        if i == rows - 1:
            ax.set_xlabel("$\zeta_{t\parallel} / \zeta_{sp}$")
        ax.set_ylabel("Fraction of bins")
        ax.set_title("$\zeta_{t\perp} = %.2f$, $|\zeta_{sp}| = %.2f$, $\\bar{n} = %i$, trials = %i" %
                     (avg_zt_y, avg_abs_zsp_x, round(expext_mean_n[i]), trials_number[i]))
        ax.legend(loc="lower right")
    plt.ion()

    fig.tight_layout()
    # fig.subplots_adjust(hspace=0.25)
    plt.show()

    plt_basename = "sim_performance"
    fig.savefig(plt_basename + ".png")
    fig.savefig(plt_basename + ".pdf")

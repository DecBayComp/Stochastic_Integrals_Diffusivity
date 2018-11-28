"""
Plot statistical data together with theoretical predictions.
Since the horizontal axis features the abs. value of zeta_sp, the regions are duplicated for the opposite signs.
Maybe I should simplify the presentation, recalculate and show the signed values? This would reduce the number of peaks and simplify presentation.
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from constants import *
from set_figure_size import set_figure_size


def plot_for_article(ksis_unique, avg_data, data_lgB, expext_mean_n, trials_number):

    # Constants
    label_location = [0.025, 0.9]
    font_size = 8
    rows = 2
    cols = 1
    dot_color = [0.1008,    0.4407,    0.7238]
    line_color = 'k'  # np.asarray([141,  24, 26]) / 255.0
    markersize = 2
    linewidth = 1
    alpha = 0.6
    xlims = [-2, 10]
    K_max = 10
    max_points = int(1e4)

    fig, _ = set_figure_size(num=1, rows=rows, page_width_frac=0.5, height_factor=2.1)
    _, ax_arr = plt.subplots(rows, cols, num=1, sharex=False, sharey=True)
    np.random.seed(0)

    for i, ax in enumerate(ax_arr):

        lgB_inferred = data_lgB[i].lg_B
        lgB_expected = data_lgB[i].lg_B_expected
        M = len(lgB_inferred)

        # Choose random indices if more than max and reduce the number of points
        indices = np.arange(M)
        if M > max_points:
            np.random.shuffle(indices)
            indices = indices[:max_points]

        lgB_inferred = lgB_inferred[indices]
        lgB_expected = lgB_expected[indices]

        plot_lims = [min(np.min(lgB_inferred), np.min(lgB_expected)),
                     max(np.max(lgB_inferred), np.max(lgB_expected))]

        # Make a linear data fit
        p = np.polyfit(lgB_expected, lgB_inferred, 1)
        y_fit = p[0] * np.array(plot_lims) + p[1]

        avg_zt_y = np.mean(avg_data[i].zeta_t_y_mean)
        if np.abs(avg_zt_y) < 1e-3:
            avg_zt_y = 0

        # % Plot
        # Fit
        l_fit = ax.plot(plot_lims, y_fit, '--', linewidth=linewidth,
                        color=line_color, label='fit')[0]
        # Identity
        l_idnt = ax.plot(plot_lims, plot_lims, color=line_color,
                         linewidth=linewidth, label='theory')[0]
        # Data
        ax.scatter(lgB_expected, lgB_inferred, marker='.',
                   s=markersize, alpha=alpha, c=dot_color)

        # Add a label to each plot
        str_label = chr(ord('a') + i)
        ax.text(label_location[0], label_location[1],
                str_label, transform=ax.transAxes, fontsize=font_size)

        ax.set_xlabel("expected $\log_{10}K^M$")
        ax.set_ylabel("inferred $\log_{10}\hat K^M$")
        ax.set_title('2D, $\zeta_{{t\perp}} = {ztper:.2f}$, bins = {M:d}'.format(
            ztper=avg_zt_y, M=max_points))

        ax.set_xlim(xlims)
        ax.set_ylim(xlims)
        ax.legend(handles=(l_idnt, l_fit), loc="lower right")

    plt.ion()

    fig.tight_layout()
    plt.show()

    plt_basename = "sim_performance"
    fig.savefig(plt_basename + ".png")
    fig.savefig(plt_basename + ".pdf")

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

matplotlib.use('Agg')


def plot_for_article(ksis_unique, avg_data, data_lgB, expext_mean_n, trials_number):

    # Constants
    label_location = [0.025, 0.9]
    font_size = 8
    rows = 1
    cols = 2
    dot_color = 'k'  # [0.1008,    0.4407,    0.7238]
    line_color = color_sequence[2]  # np.asarray([141,  24, 26]) / 255.0
    avg_color = color_sequence[1]
    avg_width = 2
    avg_alpha = 0.3
    markersize = 2
    linewidth = 1.5
    alpha = 0.6
    dot_alpha = 0.1 * 0
    xlims = [-2, 10]
    K_max = 10
    max_points = int(1e4)
    lgB_window = 0.5
    quantile = 0.025

    fig, _ = set_figure_size(num=1, rows=rows, page_width_frac=1, height_factor=0.3)
    _, ax_arr = plt.subplots(rows, cols, num=1, sharex=False, sharey=False)
    np.random.seed(0)

    for i, ax in enumerate(ax_arr):

        lgB_inferred = data_lgB[i].lg_B
        lgB_expected = data_lgB[i].lg_B_expected
        M = len(lgB_inferred)

        # # Choose random indices if more than max and reduce the number of points
        # indices = np.arange(M)
        # if M > max_points:
        #     np.random.shuffle(indices)
        #     indices = indices[:max_points]
        #
        # lgB_inferred = lgB_inferred[indices]
        # lgB_expected = lgB_expected[indices]

        lgB_min, lgB_max = [min(np.min(lgB_inferred), np.min(lgB_expected)),
                            max(np.max(lgB_inferred), np.max(lgB_expected))]
        plot_lims = [lgB_min, lgB_max]

        # Make a linear data fit
        p = np.polyfit(lgB_expected, lgB_inferred, 1)
        y_fit = p[0] * np.array(plot_lims) + p[1]

        avg_zt_y = np.mean(avg_data[i].zeta_t_y_mean)
        if np.abs(avg_zt_y) < 1e-3:
            avg_zt_y = 0

        # Get a sliding window average
        lgB_range = lgB_max - lgB_min
        steps = int(np.ceil(lgB_range / lgB_window))
        avgs = np.zeros(steps) * np.nan
        xs = np.zeros(steps) * np.nan
        CIs = np.zeros((2, steps)) * np.nan
        # lgB_inferred_sorted = np.sort(lgB_inferred)
        # lgB_expected_sorted = np.sort(lgB_expected)
        for step in range(steps):
            start = lgB_min + (step - 1) * lgB_window
            end = start + lgB_window
            xs[step] = start + lgB_window / 2
            win_data = lgB_inferred[(lgB_expected >= start) & (lgB_expected < end)]
            avgs[step] = np.mean(win_data)
            if len(win_data) > 0:
                CIs[:, step] = np.quantile(win_data, [quantile, 1 - quantile])
            # % Plot
            # Fit
        # print(CIs)

        # %% Plot
        # Data
        if dot_alpha > 0:
            ax.scatter(lgB_expected, lgB_inferred, marker='.',
                       s=markersize, alpha=dot_alpha, c=dot_color)
        # Identity
        l_idnt = ax.plot(plot_lims, plot_lims, '--', color=line_color,
                         linewidth=linewidth, label='theory')[0]
        # Avg
        l_fit = ax.plot(xs, avgs, '-', linewidth=avg_width,
                        color=avg_color, label='avg')[0]
        ax.fill_between(xs, CIs[1, :], CIs[0, :], alpha=avg_alpha, color=avg_color)

        # # Add a label to each plot
        # str_label = chr(ord('a') + i)
        # ax.text(label_location[0], label_location[1],
        #         str_label, transform=ax.transAxes, fontsize=font_size)

        ax.set_xlabel("$\log_{10}K^M$")
        ax.set_ylabel("$\log_{10}\hat K^M$")
        # ax.set_title('2D, $\zeta_{{t\perp}} = {ztper:.2f}$, trials = {tr:d}, $\\bar n = {n:d}, \Delta = {lgB_window:.2f}$'.format(
        #     ztper=avg_zt_y, M=max_points, lgB_window=lgB_window, tr=int(trials_number[i]), n=int(expext_mean_n[i])))
        ax.set_title('$\zeta_{{t\perp}} = {ztper:.2f}$'.format(
            ztper=avg_zt_y))

        ax.set_xlim(xlims)
        ax.set_ylim(xlims)
        ax.legend(handles=(l_fit, l_idnt), loc="lower right")

    plt.ion()

    fig.tight_layout()
    plt.show()

    plt_basename = "sim_performance"
    fig.savefig(plt_basename + ".png")
    fig.savefig(plt_basename + ".pdf", transparent=False)
    fig.savefig(plt_basename + ".eps", transparent=False)

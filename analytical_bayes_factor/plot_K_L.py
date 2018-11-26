

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from constants import *
from set_figure_size import set_figure_size


def plot_K_L(zeta_sps, zeta_t_roots, ns, us, dim, zeta_t_perp, height_factor=1.0, labels=True, num=1):
    """
    Make the parametric plot of the local Bayes factor
    """
    # Constants
    alpha = 0.18
    font_size = 8
    label_location = [0.05, 0.825]
    # pagewidth_in = 6.85
    # dpi = 120
    # rows = 2
    cols = len(ns)

    lambs_count = np.size(zeta_t_roots, 2)
    ns_count = len(ns)
    us_count = len(us)
    if dim == 1:
        rows = us_count
        ind_ys = range(us_count)
        ys_count = us_count

        def y_val(i): return us[i]
    elif dim == 2:
        ztpers_count = len(zeta_t_perp)
        rows = ztpers_count
        ind_ys = range(ztpers_count)
        ys_count = ztpers_count

        def y_val(i): return zeta_t_perp[i]
    dim_str = "%iD" % (dim)

    fig = set_figure_size(
        num=num, rows=rows, page_width_frac=0.5, height_factor=height_factor)
    fig, axarr = plt.subplots(rows, cols, num=num, sharey=True, sharex=True)
    # print(axarr)
    if ns_count == 1 and ys_count == 1:
        one_figure = True
        str_location = '_main'
        axarr = np.asarray([[axarr]])
    else:
        one_figure = False
        str_location = '_appendix'

    count = 0
    for ind_y in ind_ys:

        # y title
        if dim == 1:
            y_str = '$\zeta_\mathrm{t}$'
            info_str = '$u='
        elif dim == 2:
            y_str = '$\zeta_{\mathrm{t}\parallel}$'
            info_str = '$|\zeta_{t\perp}|='
        info_str += '%.1f$' % (y_val(ind_y))

        if not one_figure:
            y_str += ' (' + info_str + ')'
        axarr[ind_y, 0].set_ylabel(y_str)

        for n_ind in range(ns_count):
            ax = axarr[ind_y, n_ind]
            n = ns[n_ind]

            # Plot fixed-lambda
            for lamb_ind in range(lambs_count - 1):
                ax.fill_between(zeta_sps, zeta_t_roots[ind_y, n_ind, lamb_ind, 0, :],
                                zeta_t_roots[ind_y, n_ind, lamb_ind, 1, :],
                                color=color_sequence[lamb_ind], alpha=alpha)

            lamb_ind = lambs_count - 1
            # Positive zeta_t branch
            branch_ind = 0
            y_axis = zeta_t_roots[ind_y, n_ind, lamb_ind, branch_ind, :]
            ax.plot(zeta_sps, y_axis, color=color_sequence[lamb_ind])

            # Negative zeta_t branch
            branch_ind = 1
            y_axis = zeta_t_roots[ind_y, n_ind, lamb_ind, branch_ind, :]
            ax.plot(zeta_sps, y_axis, color=color_sequence[lamb_ind])

            # Adjust
            # plt.axis([zeta_t_points[0] * 1.05, zeta_t_points[-1] * 1.05, -2, 2])
            if ind_y == ys_count - 1:
                ax.set_xlabel('$\zeta_\mathrm{sp}$')

            if ind_y == 0:
                str_title = "$n = %i$" % (n)
                if one_figure:
                    str_title = dim_str + ', ' + str_title + ', ' + info_str
                ax.set_title(str_title)

            # Add a label to each plot
            if labels:
                str_label = chr(ord('a') + count)
                ax.text(label_location[0], label_location[1],
                        str_label, transform=ax.transAxes, fontsize=font_size)
                count += 1

    # Title
    if not one_figure:
        # dim_str += " ($\zeta_{t\perp}=%.2f$)" % (zeta_t_perp)
        fig.suptitle(dim_str, fontsize=font_size + 2)

    # plt.axis('square')
    fig.tight_layout()
    fig.subplots_adjust(top=0.89)

    # plt.subplots_adjust(wspace=0, hspace=0, left = 0.0, right = 1.0)
    fig.show()

    fig_basename = "results_%iD" % dim + str_location
    fig.savefig(fig_basename + ".pdf")
    fig.savefig(fig_basename + ".png")
    # print(plt.rcParams)

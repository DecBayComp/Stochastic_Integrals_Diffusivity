

from constants import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from set_figure_size import set_figure_size


def plot_K_L_2D(zeta_sps, zeta_t_roots, ns, us, dim, ztpers):
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
    ztpers_count = len(ztpers)
    rows = ztpers_count

    fig = set_figure_size(
        num=1, rows=rows, page_width_frac=0.5, height_factor=1.0)
    fig, axarr = plt.subplots(rows, cols, num=1, sharey=True, sharex=True)

    count = 0
    for ztpers_ind in range(ztpers_count):
        ztper = ztpers[ztpers_ind]

        # y title
        if dim == 1:
            y_str = '$\zeta_\mathrm{t}$ $(u='
        elif dim == 2:
            y_str = '$\zeta_{\mathrm{t}\parallel}$ $(|\zeta_{t\perp}|='
        y_str += '%.1f)$' % (ztper)
        axarr[ztpers_ind, 0].set_ylabel(y_str)

        for n_ind in range(ns_count):
            ax = axarr[ztpers_ind, n_ind]
            n = ns[n_ind]

            # Plot fixed-lambda
            for lamb_ind in range(lambs_count - 1):
                ax.fill_between(zeta_sps, zeta_t_roots[ztpers_ind, n_ind, lamb_ind, 0, :],
                                zeta_t_roots[ztpers_ind,
                                             n_ind, lamb_ind, 1, :],
                                color=color_sequence[lamb_ind], alpha=alpha)

            lamb_ind = lambs_count - 1
            # Positive zeta_t branch
            branch_ind = 0
            y_axis = zeta_t_roots[ztpers_ind, n_ind, lamb_ind, branch_ind, :]
            # print(y_axis)
            ax.plot(zeta_sps, y_axis, color=color_sequence[lamb_ind])

            # Negative zeta_t branch
            branch_ind = 1
            y_axis = zeta_t_roots[ztpers_ind, n_ind, lamb_ind, branch_ind, :]
            ax.plot(zeta_sps, y_axis, color=color_sequence[lamb_ind])

            # Adjust
            # plt.axis([zeta_t_points[0] * 1.05, zeta_t_points[-1] * 1.05, -2, 2])
            if ztpers_ind == ztpers_count - 1:
                ax.set_xlabel('$\zeta_\mathrm{sp}$')

            if ztpers_ind == 0:
                ax.set_title("n = %i" % (n))

            # Add a label to each plot
            str_label = chr(ord('a') + count)
            ax.text(label_location[0], label_location[1],
                    str_label, transform=ax.transAxes, fontsize=font_size)
            count += 1

    # Title
    dim_str = "%iD" % (dim)
    # if dim == 2:
    #     dim_str += " ($\zeta_{t\perp}=%.2f$)" % (zeta_t_perp)
    # # fig.text(0.01, 0.97, dim_str, fontsize = font_size + 2, weight = 'bold')
    fig.suptitle(dim_str, fontsize=font_size + 2)

    # plt.axis('square')
    fig.tight_layout()
    fig.subplots_adjust(top=0.89)

    # plt.subplots_adjust(wspace=0, hspace=0, left = 0.0, right = 1.0)
    fig.show()

    fig_basename = "results_%iD" % dim
    fig.savefig(fig_basename + ".pdf")
    fig.savefig(fig_basename + ".png")
    # print(plt.rcParams)

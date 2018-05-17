

from constants import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from set_figure_size import set_figure_size


def plot_K_L(zeta_sps, zeta_t_roots, ns, us, dim, zeta_t_perp):
    """
    Make the parametric plot of the local Bayes factor
    """
    # Constants
    alpha = 0.18
    font_size = 8
    # pagewidth_in = 6.85
    # dpi = 120
    # rows = 2
    cols = len(ns)

    lambs_count = np.size(zeta_t_roots, 2)
    ns_count = len(ns)
    us_count = len(us)
    rows = us_count

    fig = set_figure_size(
        num=1, rows=rows, page_width_frac=0.5, height_factor=1.0)
    fig, axarr = plt.subplots(rows, cols, num=1, sharey=True, sharex=True)

    # Print dimensions
    dim_str = "%iD" % (dim)
    if dim == 2:
        dim_str += " ($\zeta_{t\perp}=%.2f$)" % (zeta_t_perp)
    # fig.text(0.01, 0.97, dim_str, fontsize = font_size + 2, weight = 'bold')
    fig.suptitle(dim_str, fontsize=font_size + 2)

    for u_ind in range(us_count):
        u = us[u_ind]

        # y title
        if dim == 1:
            y_str = '$\zeta_\mathrm{t}$ $(u='
        elif dim == 2:
            y_str = '$\zeta_{\mathrm{t}\parallel}$ $(u_\perp='
        y_str += '%.1f)$' % (u)
        axarr[u_ind, 0].set_ylabel(y_str)

        for n_ind in range(ns_count):
            ax = axarr[u_ind, n_ind]
            n = ns[n_ind]

            # Plot fixed-lambda
            for lamb_ind in range(lambs_count - 1):
                ax.fill_between(zeta_sps, zeta_t_roots[u_ind, n_ind, lamb_ind, 0, :],
                                zeta_t_roots[u_ind, n_ind, lamb_ind, 1, :],
                                color=color_sequence[lamb_ind], alpha=alpha)

            lamb_ind = lambs_count - 1
            # Positive zeta_t branch
            branch_ind = 0
            y_axis = zeta_t_roots[u_ind, n_ind, lamb_ind, branch_ind, :]
            ax.plot(zeta_sps, y_axis, color=color_sequence[lamb_ind])

            # Negative zeta_t branch
            branch_ind = 1
            y_axis = zeta_t_roots[u_ind, n_ind, lamb_ind, branch_ind, :]
            ax.plot(zeta_sps, y_axis, color=color_sequence[lamb_ind])

            # Adjust
            # plt.axis([zeta_t_points[0] * 1.05, zeta_t_points[-1] * 1.05, -2, 2])
            if u_ind == us_count - 1:
                ax.set_xlabel('$\zeta_\mathrm{sp}$')

            if u_ind == 0:
                ax.set_title("n = %i" % (n))

    # plt.axis('square')
    fig.tight_layout()
    fig.subplots_adjust(top=0.89)

    # plt.subplots_adjust(wspace=0, hspace=0, left = 0.0, right = 1.0)
    fig.show()

    fig_basename = "results_%iD" % dim
    fig.savefig(fig_basename + ".pdf")
    fig.savefig(fig_basename + ".png")
    # print(plt.rcParams)



from get_v_func import get_v_func
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from set_figure_size import set_figure_size


def plot_zeta_t_perp(n_pi, us):

    # Constants
    # pagewidth_in = 6.85
    # dpi = 120
    font_size = 8
    n_min = 10
    n_max = 100
    dim = 2  # this calculation only makes sense in 2D
    B_threshold = 10
    style_sequence = ['--', '-', ':']

    # Calculate
    ns = np.asarray(range(n_min, n_max + 1))
    ns_count = len(ns)
    etas = np.sqrt(n_pi / (ns + n_pi))
    ps = dim * (ns + n_pi - 1.0) / 2.0 - 1.0

    us_count = len(us)
    zeta_t_perp = np.zeros((us_count, ns_count))

    for u_ind in range(us_count):
        u = us[u_ind]
        v = get_v_func(ns=ns, n_pi=n_pi, us=u, dim=dim, zeta_t_perp=0)

        upstairs = 1 - np.power(B_threshold / etas, 1 / ps)
        downstairs = np.multiply(
            np.power(B_threshold / etas, 1 / ps), etas ** 2) - 1
        # print(np.stack((upstairs, downstairs), 1))
        res = np.divide(np.multiply(v(0), upstairs), downstairs)
        # print(res)
        res[res < 0] = np.nan
        zeta_t_perp[u_ind, :] = np.sqrt(res)
    # print(zeta_t_perp)

    ffig = set_figure_size(num=2, rows=1, page_width_frac=0.5)
    fig, ax = plt.subplots(num=2)
    # fig.set_dpi(dpi)

    # print(ax)
    str_legend = []
    for u_ind in range(us_count):
        str_label = "u = %.1f" % us[u_ind]
        ax.plot(ns, zeta_t_perp[u_ind, :],
                linestyle=style_sequence[u_ind], label=str_label)

    # Adjust
    # ax.set_xlim(0)
    ax.set_ylim(0)
    # plt.ticklabel_format(style = "plain")
    ax.set_xlabel("n")
    ax.set_ylabel("$|\zeta_{t\perp}|$")
    ax.legend(fontsize=font_size, loc=0)

    # fig, axarr = plt.subplots(us_count, cols, num = 2, sharey = True, sharex = False)

    fig.tight_layout()
    fig.show()

try:
    has_run
except NameError:
    %matplotlib
    %load_ext autoreload
    %autoreload 2

    has_run = 1
else:
    print("Graphic interface NOT re-initialized")

import matplotlib.pyplot as plt
import numpy as np

from set_figure_size import set_figure_size

# %%
height_factor = 2 / 3
D_min = 0.01  # in um^2/s
D_ratio = 2.0
L = 1  # um
xlims = [0, L]
steps = 3
line_width = 2.5
markersize = 5
color = 'rebeccapurple'


def plot_D_profile(D_inferred_data):

    # %%
    D_max = D_min * D_ratio
    x_mesh = np.linspace(0, L, num=steps)

    def D_func(x):
        return D_max - (D_max - D_min) * np.abs(x - 0.5) * 2.0

    D_profile = D_func(x_mesh)
    D_errorbars = np.ones((2, len(D_inferred_data))) * np.nan

    # Calculation code for d errorbars (relative width)
    sorted_data = D_inferred_data.sort_values('x_center')
    x_data = sorted_data.x_center.values
    D_data = sorted_data.D.values
    D_errorbars[0, :] = D_data - sorted_data.D_CI_lower
    D_errorbars[1, :] = sorted_data.D_CI_upper - D_data

    # %%
    fig, _ = set_figure_size(num=1, rows=1, page_width_frac=0.5, height_factor=height_factor)
    _, ax = plt.subplots(1, 1, num=1, clear=True)
    ax.plot(x_mesh, D_profile, lw=line_width)

    # Experimental data
    # print(len(x_mesh), x_mesh, D_inferred_data.x_center.values,
    #       D_inferred_data.D.values, D_errorbars)

    ax.errorbar(x_data, D_data,
                yerr=D_errorbars, c=color, capsize=2, zorder=3, ms=markersize, fmt='x')

    # plt.title('{n} bins'.format(n=len(D_inferred_data.D)))

    ylims = [0, D_max * 1.025]
    plt.xlim(xlims)
    plt.ylim(0)
    plt.xticks([0, 0.5, 1])
    plt.yticks([0, 0.01, 0.02])
    plt.xlabel('$x, \mu \mathrm{m}$')
    plt.ylabel('$D, \mu \mathrm{m^2/s}$')

    fig.tight_layout()
    plt.show()

    figname = 'D_profile'
    fig.savefig(figname + '.pdf')
    fig.savefig(figname + '.png')
    fig.savefig(figname + '.eps')

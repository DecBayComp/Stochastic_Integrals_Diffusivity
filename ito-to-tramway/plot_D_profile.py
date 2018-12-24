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
rows = 1
cols = 1
height_factor = 2
D_min = 0.01  # in um^2/s
D_ratio = 2.0
L = 1  # um
xlims = [0, L]
steps = 3
line_width = 2.5


# %%
D_max = D_min * D_ratio
x_mesh = np.linspace(0, L, num=steps)


def D_func(x):
    return D_max - (D_max - D_min) * np.abs(x - 0.5) * 2.0


D_profile = D_func(x_mesh)

# %%
fig, _ = set_figure_size(num=1, rows=rows, page_width_frac=0.5, height_factor=height_factor)
_, ax_arr = plt.subplots(rows, cols, num=1, sharex=False, sharey=True)
plt.plot(x_mesh, D_profile, lw=line_width)

ylims = [0, D_max * 1.025]
plt.xlim(xlims)
plt.ylim(ylims)
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

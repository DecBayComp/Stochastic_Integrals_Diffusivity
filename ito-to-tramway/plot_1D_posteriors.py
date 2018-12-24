"""
Plot 4 posteriors for the given bin Parameters
"""

try:
    has_run
except NameError:
    %matplotlib
    %load_ext autoreload
    %autoreload 2

    has_run = 1
else:
    print("Graphic interface NOT re-initialized")

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from constants import color_sequence
from set_figure_size import set_figure_size
from stopwatch import stopwatch
from tramway.inference.bayes_factors import calculate_bayes_factors
from tramway.inference.bayes_factors.calculate_posteriors import (calculate_one_1D_posterior_in_2D,
                                                                  calculate_one_1D_prior_in_2D)

# %% Estimate the expected zetas
D_min = 0.01  # in um^2/s
D_ratio = 2.0
D_max = D_min * D_ratio
L = 1  # um
steps = 500
ksi = 0.5
zty = 0.1 * 0
u = 1
D_grad_abs = (D_max - D_min) / L * 2
V = 1  # 2 * dim * D * dt
std = np.sqrt(V)
V_pi = u * V
loc_error = 0.1**2 * 0
lamb_list = [0, 0.5, 1, 'marg']
V = 1
V_pi = 1
axis = 'x'
rows = 1
cols = 1
height_factor = 2
bl_appendix = True

if not bl_appendix:
    zeta_sp = [0.1, 0]
    zeta_t = [0, 0]
    zeta_a_lims = [-0.5, 0.5]
    n = 500
    figname = '4_posteriors'
else:
    zeta_sp = [0.5, 0]
    zeta_t = [0, 0]
    zeta_a_lims = [-1, 1]
    n = 100
    figname = 'appendix_4_posteriors_skewed'


def D_func(x):
    return D_max - (D_max - D_min) * np.abs(x - 0.5) * 2.0


# %% Calculate posteriors
zax_mesh = np.linspace(zeta_a_lims[0], zeta_a_lims[1], num=steps)
p_zax = {}  # = np.zeros(steps) * np.nan
with stopwatch():
    for lamb in lamb_list:
        p_zax[lamb] = np.zeros(steps) * np.nan
        for i, zax in enumerate(tqdm(zax_mesh)):
            p_zax[lamb][i] = calculate_one_1D_posterior_in_2D(
                zeta_a=zax, zeta_t=zeta_t, zeta_sp=zeta_sp, n=n, V=V, V_pi=V_pi, loc_error=loc_error, lamb=lamb, axis=axis)

# Calculate the prior
lamb = 'prior'
p_zax[lamb] = np.zeros(steps) * np.nan
for i, zax in enumerate(tqdm(zax_mesh)):
    p_zax[lamb][i] = calculate_one_1D_prior_in_2D(
        zeta_a=zax, V_pi=V_pi, loc_error=loc_error)


# lg_BM = calculate_bayes_factors(zeta_ts=np.array([zeta_t]), zeta_sps=np.array([zeta_sp]), ns=[n],
#                                 Vs=[V], Vs_pi=[V_pi], loc_error=loc_error, dim=2, verbose=False)[0][0]

# Detect the ROI
threshold = 1e-2
min = zax_mesh[np.argmax(p_zax[lamb_list[0]])]
max = min
for key in list(p_zax.keys())[:-1]:
    z1 = zax_mesh[p_zax[key] >= threshold]
    min = np.min([np.min(z1), min])
    max = np.max([np.max(z1), max])
zax_lims = [min, max]


# %% Plot
fig, _ = set_figure_size(num=2, rows=rows, page_width_frac=0.5, height_factor=height_factor)
_, ax_arr = plt.subplots(rows, cols, num=2, sharex=False, sharey=True)
# Prior
line_width = 1
line_type = ':'
plt.plot(zax_mesh, p_zax['prior'], line_type, color='k', linewidth=line_width)
# Fixed lambda
line_width = 1
line_type = '--'
plt.plot(zax_mesh, p_zax[0], line_type, color=color_sequence[0], linewidth=line_width)
plt.plot(zax_mesh, p_zax[0.5], line_type, color=color_sequence[1], linewidth=line_width)
plt.plot(zax_mesh, p_zax[1], line_type, color=color_sequence[2], linewidth=line_width)
# Marginalized
line_width = 2
line_type = '-'
plt.plot(zax_mesh, p_zax['marg'], line_type, color=color_sequence[3], linewidth=line_width)

# Adjust
plt.xlabel('$\zeta_{a\parallel}$')
plt.ylabel('$p(\zeta_{a\parallel} \mid T)$')
plt.xlim(zax_lims)

str_title = '2D, $\zeta_{{sp}} = {zsp}$'.format(zsp=zeta_sp[0])
str_title += ', $\zeta_{{t\parallel}}/\zeta_{{sp}} = {ksi}$'.format(ksi=zeta_t[0] / zeta_sp[0])
str_title += ', $\zeta_{{t\perp}} = {zty}$'.format(zty=zeta_t[1])
str_title += ', $n = {n}$'.format(n=n)
# str_title += ', $, \log_{{10}} B^M={lg_BM:.2f}$'.format(lg_BM=lg_BM)
plt.title(str_title)

fig.tight_layout()
plt.show()

fig.savefig(figname + '.pdf')
fig.savefig(figname + '.png')
fig.savefig(figname + '.eps')

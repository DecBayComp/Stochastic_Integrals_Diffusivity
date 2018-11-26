"""
Plot analytical Bayes factor behavior for the main text
"""

try:
    has_run
except NameError:
    %matplotlib
    %load_ext autoreload
    %autoreload 2
    has_run = True


import copy
import logging

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D, axes3d
from tqdm import tqdm, trange

from calculate_bayes_factor_limits import calculate_bayes_factor_limits
from convenience_functions import n_pi_func
from set_figure_size import set_figure_size
from tramway.inference.bayes_factors import calculate_bayes_factors

# %% Main calculations

# Constants
zeta_sp_abs_lim = 1.00  # 0.02
N = 100
lgB_levels = [-1, 1]  # at least strong evidence against H0
zeta_t_perp_2D = 0.1 * 0
n = 500
u = 1.0
dim = 2
eta = np.sqrt(n_pi_func(dim) / (n + n_pi_func(dim)))

# Meshes
zsps = np.linspace(-zeta_sp_abs_lim, zeta_sp_abs_lim, num=N)
zts = zsps

# %% Calculate
lg_Bs = np.zeros((N, N)) * np.nan
for ix in trange(N):
    for iy in range(N):
        zsp_vec = np.array([[zsps[ix], 0]])
        zt_vec = np.array([[zts[iy], zeta_t_perp_2D]])
        lg_Bs[iy, ix] = calculate_bayes_factors(zeta_ts=zt_vec, zeta_sps=zsp_vec, ns=[n],
                                                Vs=[1], Vs_pi=[u], loc_error=0, dim=dim, verbose=False)[0][0]


# %% Plot
# Do not show the colormap for the outside regions
cut_lim = -1 + 0.01
lg_Bs_cut = copy.deepcopy(lg_Bs)
lg_Bs_cut[lg_Bs_cut > cut_lim] = cut_lim * np.nan

fig = set_figure_size(num=1, rows=1, page_width_frac=0.5, height_factor=2.4)
ax = fig.subplots()

# Colormesh + a nice colorbar
c = ax.pcolormesh(zsps, zts, lg_Bs_cut, cmap='inferno', clim=[-5, 5])

gca_bkp = plt.gca()
divider = make_axes_locatable(fig.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
fig.colorbar(c, cax=cax)
plt.sca(gca_bkp)
colorbar_legend = '$\log_{10}K$'
cax.set_ylabel(colorbar_legend, rotation=90)

# Contours for the main B thresholds
color = 'k'
c_major = ax.contour(zsps, zts, lg_Bs, levels=lgB_levels, colors=color, linewidths=1)

# Outside contours
start = 10
step = 40
linewidth = 0.5
alpha = 0.3
label_coords = [(-0.8434066225109812, 0.22193196803229043),
                (-0.8510420449603688, 0.5131273273250199),
                (-0.8510415968800961, 0.7225019254197504),
                (-0.8663145118087212, 0.9133730114026575),
                (0.8748413979235168, -0.22191754989132562),
                (0.8595704914313069, -0.5131270515852441),
                (0.8595705741642174, -0.7225011437390221),
                (0.8519360978457273, -0.9133748850667568)]
lgB_minor_levels_out = np.arange(start, np.ceil(np.max(lg_Bs)) + step, step=step)
c_minor_out = ax.contour(zsps, zts, lg_Bs, levels=lgB_minor_levels_out,
                         colors=color, alpha=alpha, linestyles='solid', linewidths=linewidth)
out_labels = ax.clabel(c_minor_out, fmt='%i', manual=label_coords)
label_coords = [lab.get_position() for lab in out_labels]


# # Add minor B levels
# step = -0.2
# alpha = 0.6
# label_coords_in = [(-0.5743244447928386, -0.29449109602137313),
#                    (0.5750595136833536, 0.26666617912833934)]
# lgB_minor_levels_in = np.arange(-1 + step, np.floor(np.min(lg_Bs)) + step, step=step)
# lgB_minor_levels_in = lgB_minor_levels_in[::-1]
# c_minor_in = ax.contour(zsps, zts, lg_Bs, levels=lgB_minor_levels,
#                         colors='w', alpha=alpha, linestyles='solid', linewidths=1)
# in_labels = ax.clabel(c_minor_in, manual=label_coords_in, fmt='%1.1f', rightside_up=True)
# label_coords_in = [lab.get_position() for lab in in_labels]
# ax.clabel(c_minor, manual=0, fmt='%1.1f', rightside_up=True)
# # Their labels
# lgB_minor_labels_lims = [-1.45, -1.25]
# lgB_minor_labels = np.array([i for i in c_minor.levels if (
#     i >= lgB_minor_labels_lims[0]) & (i <= lgB_minor_labels_lims[1])])

# Labels and title
plt.xlabel('$\zeta_\mathrm{sp}$')
plt.ylabel('$\zeta_{\mathrm{t}\parallel}$')
str_title = '2D, $n={n}$, $u={u}$, $\zeta_{{t\perp}}={ztper}$'.format(
    n=n, u=u, ztper=zeta_t_perp_2D)
plt.title(str_title)
fig.tight_layout()

# Save
name = 'results_2D_main'
plt.savefig(name + '.pdf')
plt.savefig(name + '.png')


# # %% 3D figure
# lg_Bs_cut = copy.deepcopy(lg_Bs)
# cut_lim = 1
# lg_Bs_cut[lg_Bs_cut > cut_lim] = cut_lim
# # colors = ['yellow']
# fig = plt.figure(num=2, clear=True)
# ax = fig.gca(projection='3d')
# # c = plt.imshow(lg_Bs, cmap='inferno', origin='lower', interpolation='bilinear')
# X, Y = np.meshgrid(zsps, zts)
# c = ax.plot_surface(X, Y, lg_Bs_cut, cmap=cm.inferno)
#
# # # Axes
# # linewidth = 0.5
# # alpha = 0.3
# # plt.plot(plt.xlim(), [0, 0], 'w', linewidth=linewidth, alpha=alpha)
# # plt.plot([0, 0], plt.ylim(), 'w', linewidth=linewidth, alpha=alpha)
#
# # # Add contours for B thresholds
# # c2 = ax.contour(zsps, zts, lg_Bs, levels=lgB_levels, colors=colors)
#
#
# # plt.colorbar(c)
# plt.xlabel('$\zeta_\mathrm{sp}$')
# plt.ylabel('$\zeta_{\mathrm{t}\parallel}$')
#
#
# # %%
# zsp_vec = np.array([[0.5, 0]])
# zt_vec = np.array([[0.25, 0]])
# calculate_bayes_factors(zeta_ts=zt_vec, zeta_sps=zsp_vec, ns=n, Vs=[1], Vs_pi=[u],
#                         loc_error=0, dim=dim, verbose=False)
#

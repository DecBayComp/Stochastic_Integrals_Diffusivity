"""
This file regroups different specific plotting functions
"""

import logging

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

linewidth = 0.1

pagewidth_in = 6.85
font_size = 8
dpi = 100
page_width_frac = 1 / 3
color = 'k'
# alpha_mesh = 0.15

fig_name = 'bead_lattice'


def plot_beads(beads, figure=None, xlims=None, ylims=None, ticks=True, letter_label=False, colorbar=False, pdf=True, **kwargs):
    """
    If a figure is given, only the beads are plotted.
    [to-do] Otherwise, the plot is adjusted.
    """
    if not figure:
        # Set default figure font size and LaTeX usage
        matplotlib.rcParams.update({'font.size': font_size})
        fig, ax = plt.subplots(num=1, clear=True)
        ax.set_aspect('equal')

    patches = []
    for i in range(beads.shape[0]):
        bead = beads[i, :]
        patches.append(plt.Circle(bead[:2], bead[2]))
        #  facecolor='r', edgecolor=None

    ax = plt.gca()
    ax.add_collection(matplotlib.collections.PatchCollection(
        patches, facecolors=color, edgecolors=None))

    if not figure:
        figsize = np.asarray([3.0, 1.0]) * page_width_frac * pagewidth_in  # in inches
        figsize = tuple(figsize)
        # Manual figure adjustments
        fig = plt.gcf()
        fig.set_dpi(dpi)
        fig.set_figwidth(figsize[0])
        fig.set_figheight(figsize[1])

        if xlims:
            plt.xlim(xlims)
        if ylims:
            plt.ylim(ylims)

        # Remove ticks
        if not ticks:
            plt.tick_params(
                axis='x',           # changes apply to the x-axis
                which='both',       # both major and minor ticks are affected
                bottom=False,       # ticks along the bottom edge are off
                top=False,          # ticks along the top edge are off
                labelbottom=False)  # labels along the bottom edge are off
            plt.tick_params(
                axis='y',           # changes apply to the y-axis
                which='both',       # both major and minor ticks are affected
                left=False,         # ticks along the bottom edge are off
                right=False,        # ticks along the top edge are off
                labelleft=False)    # labels along the bottom edge are off
        else:
            plt.xlabel('$x, \mu\mathrm{m}$')
            plt.ylabel('$y, \mu\mathrm{m}$')

        # add label
        if letter_label:
            label_location = [0.025, 1.03]
            # str_label = chr(ord('a') + plot_me.count)
            # plot_me.count += 1
            ax = fig.gca()
            ax.text(label_location[0], label_location[1],
                    letter_label, transform=ax.transAxes, fontsize=font_size)

        # Colorbar legend
        if colorbar and colorbar_legend:
            fig.axes[1].set_ylabel(colorbar_legend, rotation=90)

        if pdf:
            try:
                fig.savefig(fig_name + '.pdf', bbox_inches='tight', pad_inches=0)
            except:
                logging.warn('Unable to save figure. The file might be open')

        # Save zoomed .png version
        factor = 5
        fig.set_figwidth(figsize[0] * factor)
        fig.set_figheight(figsize[1] * factor)
        try:
            fig.savefig(fig_name + '.png', pad_inches=0, bbox_inches='tight')  # )
        except:
            logging.warn('Unable to save figure. The file might be open')

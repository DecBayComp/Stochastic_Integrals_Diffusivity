

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def set_figure_size(num, rows, page_width_frac, height_factor=1.0):
    pagewidth_in = 6.85
    font_size = 8
    dpi = 100

    figsize = np.asarray([1.0, rows *
                          height_factor]) * page_width_frac * pagewidth_in  # in inches

    # Set default font size
    matplotlib.rcParams.update({'font.size': font_size})

    # Enable LaTeX and set font to Helvetica
    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble'] = [
        r'\usepackage{tgheros}',    # helvetica font
        r'\usepackage{sansmath}',   # math-font matching  helvetica
        r'\sansmath'                # actually tell tex to use it!
        r'\usepackage{siunitx}',    # micro symbols
        r'\sisetup{detect-all}',    # force siunitx to use the fonts
    ]
    # Enforce TrueType fonts for easier editing later on
    # matplotlib.rcParams['pdf.fonttype'] = 42
    # matplotlib.rcParams['ps.fonttype'] = 42

    # Create and return figure handle
    fig = plt.figure(num)
    fig.clf()

    # Set figure size and dpi
    fig.set_dpi(dpi)
    fig.set_figwidth(figsize[0])
    fig.set_figheight(figsize[1])

    # fig.tight_layout()

    # Return figure handle
    return (fig, figsize)

import matplotlib
import numpy as np

colors = ['r', 'g', 'b', 'c', 'y', 'm', 'k']

def color_iter():
    return iter(colors)
# A colorblind safe, print friendly, and photocopy safe color scheme:
# http://colorbrewer2.org/?type=diverging&scheme=PuOr&n=3
fig_orange = '#f1a340'
fig_gray = '#f7f7f7'
fig_purple = '#998ec3'

def set_fig_params_for_publication():
    matplotlib.rcParams['font.size'] = 6
    matplotlib.rcParams['font.sans-serif'] = 'Arial'
    matplotlib.rcParams['text.usetex'] = True
    matplotlib.rcParams['text.latex.preamble'] = [
            r'\usepackage{helvet}',
            r'\usepackage{sansmath}',
            r'\sansmath']
    matplotlib.rcParams['xtick.major.size'] = 2
    matplotlib.rcParams['ytick.major.size'] = 2
    matplotlib.rcParams['xtick.major.pad'] = 2
    matplotlib.rcParams['ytick.major.pad'] = 2

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def set_axis_ticks(ax):
    """Set ticks on left and bottom axis, ticks facing out."""
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.yaxis.set_tick_params(direction='out')
    ax.xaxis.set_tick_params(direction='out')

def format_axis(ax, label_padding=2):
    set_axis_ticks(ax)
    ax.xaxis.labelpad = label_padding
    ax.yaxis.labelpad = label_padding


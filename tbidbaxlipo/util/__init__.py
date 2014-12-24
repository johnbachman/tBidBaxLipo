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

fontsize = 6
capsize = 1.5

def set_fig_params_for_publication():
    matplotlib.rcParams['font.sans-serif'] = 'Arial'
    matplotlib.rcParams['text.usetex'] = True
    matplotlib.rcParams['text.latex.preamble'] = [
            r'\usepackage{helvet}',
            r'\usepackage{sansmath}',
            r'\sansmath',
            r'\usepackage{underscore}',]
    #matplotlib.rcParams['xtick.major.size'] = 2
    #matplotlib.rcParams['ytick.major.size'] = 2
    #matplotlib.rcParams['xtick.major.pad'] = 2
    #matplotlib.rcParams['ytick.major.pad'] = 2

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def format_axis(ax, label_padding=2, tick_padding=0, yticks_position='left'):
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position(yticks_position)
    ax.yaxis.set_tick_params(which='both', direction='out', labelsize=fontsize,
                             pad=tick_padding, length=2, width=0.5)
    ax.xaxis.set_tick_params(which='both', direction='out', labelsize=fontsize,
                             pad=tick_padding, length=2, width=0.5)
    ax.xaxis.labelpad = label_padding
    ax.yaxis.labelpad = label_padding
    ax.xaxis.label.set_size(fontsize)
    ax.yaxis.label.set_size(fontsize)

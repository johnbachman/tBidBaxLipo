import os
import sys
import csv
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
from itertools import product
from tbidbaxlipo.util import format_axis, set_fig_params_for_publication,\
                             fontsize
from tbidbaxlipo.plots.nbd_bax_analysis import site_region, color_dict
set_fig_params_for_publication()

evi_dict = {}
activators = set()
residues = set()
reps = set()
min_conf = 2
max_conf = 5
# The number of distinct conformational models tested
nconf_models = max_conf - min_conf + 1

# Set path to .csv file containing evidence values
curdir = os.path.dirname(__file__)
csv_filename = os.path.join(curdir, 'mcmc', 'pt_data1_evidence.csv')
print csv_filename
with open(csv_filename) as csv_file:
    csvreader = csv.reader(csv_file, delimiter=',')
    for row in csvreader:
        # Parse the row
        (act, res, rep, nconfs, evi, err) = row
        # Update our dict of evidence values
        nconfs = int(nconfs)
        evi_arr = evi_dict.setdefault((act, res, rep), np.zeros(nconf_models))
        evi_arr[nconfs - 2] = -float(evi)
        # Update our lists of keys
        activators.add(act)
        residues.add(res)
        reps.add(rep)

num_rows = len(activators) * len(residues) * len(reps)

def plot_3confs(plot_filename):
    min_conf_to_plot = 3
    fig = plt.figure(figsize=(1.5, 1.8), dpi=300)
    ax = fig.gca()

    # Create an array to hold all of the values
    evi_matrix = np.zeros((nconf_models, num_rows))
    # Now, iterate over all the keys
    row_index = 0
    xpositions = range(min_conf_to_plot, max_conf + 1)
    arr_start_index = min_conf_to_plot - min_conf
    for (act, res, rep) in product(activators, residues, reps):
        evi_arr = evi_dict[(act, res, rep)]
        evi_matrix[:, row_index] = evi_arr
        # Get jitter values for the xcoordinates
        xjitter = np.random.randn(len(xpositions)) * 0.07
        xvalues = xpositions + xjitter
        # Color code according to part of the protein
        color = color_dict[site_region[res]]
        # Plot
        ax.plot(xvalues, evi_arr[arr_start_index:], linestyle='', marker='.',
                markersize=2, color=color)
        row_index += 1

    ax.set_xlim(min_conf_to_plot - 0.5, max_conf + 0.5)
    ax.set_xticks(xpositions)
    ax.set_xticklabels([str(xpos) for xpos in xpositions])
    #ax.set_xlabel('Conformations in model')
    #ax.set_ylabel('-ln(Marginal likelihood)')

    sub_matrix = evi_matrix[arr_start_index:]
    ax.boxplot(sub_matrix.T, positions=xpositions, widths=0.8, sym='',
               boxprops={'color':'black'},
               meanprops={'color':'black'},
               medianprops={'color':'black'},
               whiskerprops={'color':'black', 'linestyle':'-'})

    format_axis(ax)

    # Create legend with lines to match colors of points
    n_patch = mlines.Line2D([], [], color=color_dict['nterm'],
                            label='N-term')
    bh3_patch = mlines.Line2D([], [], color=color_dict['bh3'],
                              label='BH3')
    a56_patch = mlines.Line2D([], [], color=color_dict['h56'],
                              label='$\\alpha$5-6')
    a9_patch = mlines.Line2D([], [], color=color_dict['h9'],
                              label='$\\alpha 9$')
    leg = plt.legend(handles=[n_patch, bh3_patch, a56_patch, a9_patch],
             loc='upper right', borderaxespad=0.,
             prop={'size': fontsize}, handlelength=2,
             handletextpad=0.2, labelspacing=0.5)
    leg.draw_frame(False)

    plt.subplots_adjust(left=0.16, bottom=0.12)
    plt.savefig('%s.pdf' % plot_filename)
    plt.savefig('%s.png' % plot_filename)

def plot_2confs(plot_filename):
    min_conf_to_plot = 2
    f,(ax1, ax2) = plt.subplots(2, 1, figsize=(3, 3), dpi=300)

    # Create an array to hold all of the values
    evi_matrix = np.zeros((nconf_models, num_rows))
    # Now, iterate over all the keys
    row_index = 0
    xpositions = range(min_conf_to_plot, max_conf + 1)
    arr_start_index = min_conf_to_plot - min_conf
    for (act, res, rep) in product(activators, residues, reps):
        evi_arr = evi_dict[(act, res, rep)]
        evi_matrix[:, row_index] = evi_arr
        # Get jitter values for the xcoordinates
        xjitter = np.random.randn(len(xpositions)) * 0.05
        xvalues = xpositions + xjitter
        # Color code according to part of the protein
        color = color_dict[site_region[res]]
        # Plot
        for ax in (ax1, ax2):
            ax.plot(xvalues, evi_arr[arr_start_index:], linestyle='', marker='.',
                    markersize=2, color=color)
        row_index += 1

    ax1.set_xlim(min_conf_to_plot - 0.5, max_conf + 0.5)
    ax2.set_xlim(min_conf_to_plot - 0.5, max_conf + 0.5)
    ax2.set_xticks(xpositions)
    ax2.set_xticklabels([str(xpos) for xpos in xpositions])
    ax2.set_xlabel('Conformations in model')

    ax1.set_ylim(1600, 30000)
    ax2.set_ylim(0, 1300)
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    #ax1.xaxis.tick_top()
    #ax1.tick_params(labeltop='off')
    #ax2.xaxis.tick_bottom()

    # Diagonal lines
    d = .015 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d,+d),(-d,+d), **kwargs)      # top-left diagonal
    ax1.plot((1-d,1+d),(-d,+d), **kwargs)    # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d,+d),(1-d,1+d), **kwargs)   # bottom-left diagonal
    ax2.plot((1-d,1+d),(1-d,1+d), **kwargs) # bottom-right diagonal

    sub_matrix = evi_matrix[arr_start_index:]
    ax2.boxplot(sub_matrix.T, positions=xpositions, widths=0.8, sym='',
               boxprops={'color':'black'},
               meanprops={'color':'black'},
               medianprops={'color':'black'},
               whiskerprops={'color':'black', 'linestyle':'-'})

    # Set up the yaxis for the top plot
    #ax1.yaxis.set_tick_params(which='both', direction='out', labelsize=fontsize,
    #                         pad=0, length=2, width=0.5)
    #ax1.yaxis.set_ticks_positions('left')
    #ax1.yaxis.labelpad = 2
    #ax1.yaxis.label.set_size(fontsize)
    #ax1.set_xticks([])
    format_axis(ax1)
    ax1.set_xticks([])
    format_axis(ax2)
    f.subplots_adjust(hspace=0.05)

    # Add shared yaxis label
    f.text(0.035, 0.5, '-ln(Marginal likelihood)', ha='center', va='center',
           rotation='vertical', fontsize=fontsize)

    plt.subplots_adjust(left=0.16, bottom=0.12)
    plt.savefig('%s.pdf' % plot_filename)
    plt.savefig('%s.png' % plot_filename)

if __name__ == '__main__':
    plot_2confs('test2confs')
    plot_3confs('test3confs')


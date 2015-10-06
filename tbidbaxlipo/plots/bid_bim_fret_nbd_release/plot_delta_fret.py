from tbidbaxlipo.plots import nbd_bax_analysis as nba
from matplotlib import pyplot as plt
from tbidbaxlipo.models.nbd.multiconf import Builder
import numpy as np
import cPickle
from tbidbaxlipo.util import set_fig_params_for_publication, format_axis, \
                             fontsize
from tbidbaxlipo.plots.bid_bim_fret_nbd_release.preprocess_data \
    import df_pre, nbd_residues
import matplotlib.patches as mpatches
import sys

def calc_fret_deltas(residues):
    delta_dict = {}

    # --- First, fit Bid data ---
    bid_fit_results = nba.plot_3conf_fits(df_pre, residues, 'Bid', dtype='FRET',
                                          do_plot=False)
    for fr in bid_fit_results:
        # Calculate the difference
        key = (fr.activator, fr.measurement, fr.nbd_site, fr.rep_index)
        delta_fret = np.max(fr.y) - fr.y[-1]
        delta_dict[key] = delta_fret

    # --- Now do Bim ---
    bim_fit_results = nba.plot_3conf_fits(df_pre, residues, 'Bim', dtype='FRET',
                                          do_plot=False)
    for fr in bim_fit_results:
        # Calculate the difference
        key = (fr.activator, fr.measurement, fr.nbd_site, fr.rep_index)
        delta_fret = np.max(fr.y) - fr.y[-1]
        delta_dict[key] = delta_fret

    # Cache the results
    with open('fret_deltas.pck', 'w') as f:
        cPickle.dump(delta_dict, f)

    return delta_dict

def load_fret_deltas():
    # Load data
    with open('fret_deltas_bid.pck') as f:
        bid_deltas = cPickle.load(f)
    with open('fret_deltas_bim.pck') as f:
        bim_deltas = cPickle.load(f)
    bid_deltas.update(bim_deltas)
    return bid_deltas

def get_keys(deltas):
    activators = set()
    residues = set()
    reps = set()
    for key in deltas.keys():
        activators.add(key[0])
        residues.add(key[2])
        reps.add(key[3])
    sorted_residues = sorted(residues, key=lambda x: int(x))
    return (sorted(activators), sorted_residues, sorted(reps))

def plot_means(deltas, plot_filename=None):
    # Calculate the means
    (activators, residues, reps) = get_keys(deltas)
    residues = [res for res in residues if res != '62']
    set_fig_params_for_publication()
    (fig_width, rel_left, rel_right) = nba.calc_barplot_width(len(residues))
    bar_colors = {'Bid': 'gray', 'Bim': 'black'}

    plt.figure(figsize=(fig_width, 1.5), dpi=300)
    for act_ix, act in enumerate(activators):
        for nbd_ix, res in enumerate(residues):
            vals = [deltas[(act, 'FRET', res, rep)] for rep in reps]
            #means[(act, res)] = np.mean(vals)
            #sds[(act, res)] = np.std(vals, ddof=1)

            plt.bar(range(nbd_ix*3 + act_ix, (nbd_ix*3) + 1 + act_ix),
                    np.mean(vals), width=1, color=bar_colors[act],
                    linewidth=0, ecolor='k', capsize=1.5,
                    yerr=np.std(vals, ddof=1))
    x_lbound = -1
    x_ubound = len(residues) * 3
    plt.hlines(0, x_lbound, x_ubound, color='0.85', zorder=1, linewidth=0.5)
    plt.subplots_adjust(left=rel_left, bottom=0.1, right=rel_right, top=0.94)
    ax = plt.gca()
    ax.set_xlim([x_lbound, x_ubound])
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_tick_params(which='both', labelsize=fontsize, pad=2,
            length=2, width=0.5)
    ax.yaxis.set_tick_params(which='both', direction='out',
                    labelsize=fontsize, pad=0, length=2, width=0.5)
    ax.xaxis.labelpad = 2
    ax.yaxis.labelpad = 2
    ax.set_xticks(np.arange(1, 1 + len(residues) * 3, 3))
    ax.set_xticklabels(residues)
    ax.set_ylabel('Max FRET - Endpoint FRET', fontsize=fontsize)
    # Create legend with rectangles to match colors of bars
    bid_patch = mpatches.Patch(color=bar_colors['Bid'], label='cBid')
    bim_patch = mpatches.Patch(color=bar_colors['Bim'], label='Bim')
    leg = plt.legend(handles=[bid_patch, bim_patch],
            bbox_to_anchor=(1.02, 0.5), loc='center left', borderaxespad=0.,
            prop={'size': fontsize}, handlelength=1)
    leg.draw_frame(False)

    # Save the figure
    if plot_filename:
        plt.savefig('%s.pdf' % plot_filename)
        plt.savefig('%s.png' % plot_filename, dpi=300)

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print("Usage: python plot_delta_fret.py plot_filename")
        sys.exit(1)
    plot_filename = sys.argv[1]

    residues = [res for res in nbd_residues if res != '62']
    deltas = calc_fret_deltas(residues)
    plot_means(deltas, plot_filename)



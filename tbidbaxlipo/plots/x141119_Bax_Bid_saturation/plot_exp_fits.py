import sys
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import linregress

from tbidbaxlipo.util.plate_assay import plot_all, TIME, VALUE
from tbidbaxlipo.util import fitting, colors, set_fig_params_for_publication, \
                             format_axis
from preprocess_data import bim_bh3, bid_80, bid_40, bid_20, bid_10, \
                            bid_5, bid_2, bid_0, bax_concs, bid_concs, \
                            bg_averages, bg_diff, data_matrix, data_norm
from tbidbaxlipo.plots import titration_fits as tf

def plot_bg():
    plt.figure()
    plot_all(bg_averages)
    plt.title('Measured background conditions')

    plt.figure()
    plot_all(bg_diff)
    plt.title('Calculated background conditions')

def plot_data():
    """Plots the data and various transformations of it."""
    plt.figure()
    plot_all(bim_bh3)
    plt.title('Bim BH3')

    plt.figure()
    plot_all(bid_80)
    plt.title('Bid 80nM')

    plt.figure()
    plot_all(bid_40)
    plt.title('Bid 40 nM')

    plt.figure()
    plot_all(bid_20)
    plt.title('Bid 20 nM')

    plt.figure()
    plot_all(bid_10)
    plt.title('Bid 10 nM')

    plt.figure()
    plot_all(bid_5)
    plt.title('Bid 5 nM')

    plt.figure()
    plot_all(bid_2)
    plt.title('Bid 2 nM')

    plt.figure()
    plot_all(bid_0)
    plt.title('Bid 0 nM')

def plot_mm(k_data, bax_concs, bid_concs):
    plt.figure('mm fit', figsize=(1.5, 1.5), dpi=300)

    kcat = fitting.Parameter(0.06)
    km = fitting.Parameter(250.)
    v0 = fitting.Parameter(5e-5)
    ekd = fitting.Parameter(2.)

    def plot_fit():
        for i, bid in enumerate(bid_concs):
            c = colors[i]
            bid_bound = bid / (ekd() + bid)
            if bid == 2.5:
                bid_str = '2.5'
            else:
                bid_str = str(int(bid))
            plt.plot(bax_concs,
                     ((kcat() * bid_bound) / (km() + bax_concs)) + v0(),
                     color=c, label='%s nM' % bid_str)
            bid_k = k_data[i]
            plt.plot(bax_concs, bid_k, marker='o', color=c,
                     linestyle='', markersize=3)

    def fit_func(bax_concs):
        res_list = []
        for bid in bid_concs:
           bid_bound = bid / (ekd() + bid)
           res_list.append(((kcat() * bid_bound) / (km() + bax_concs)) + v0())
        return np.hstack(res_list)

    fitting.fit(fit_func, [kcat, km, v0, ekd],
                np.hstack(k_data), np.array(bax_concs))

    plot_fit()

    plt.subplots_adjust(left=0.21, bottom=0.19)
    plt.xlabel('[Total Bax] (nM)')
    plt.ylabel(r'k (sec$^{-1} \times 10^{-5}$)')
    ax = plt.gca()
    ax.set_ylim([5e-5, 37e-5])
    ax.set_xlim([10, 2000])
    ax.set_yticks(np.linspace(5e-5, 35e-5, 7))
    ax.set_yticklabels([int(f) for f in np.linspace(5, 35, 7)])
    ax.set_xscale('log')
    format_axis(ax)

    #ax.legend(handles, loc='upper right', ncol=1,
    #            borderpad=0,
               #bbox_to_anchor=(0.7, 0.9),
               #bbox_transform=plt.gcf().transFigure,
    #           handlelength=1.5, handletextpad=0,
    #           frameon=False, prop={'size':6})

    #plt.text(2.5, 0.00033, '$K_{cat} = %.4f\ sec^{-1}$' % kcat())
    #plt.text(2.5, 0.00030, '$K_m = %.2f\ nM$' % km() )
    #plt.text(2.5, 0.00027, '$V_0 = %f\ sec^{-1}$' % v0())
    #plt.text(2.5, 0.00024, '$Bid/Lipo\ K_D = %.2f\ nM$' % ekd())

def calc_exp_fits(data_norm):
    """Fit each of the timecourses in the data matrix to an exponential func.

    Parameters
    ----------
    data_norm : np.array
        Four-dimensional numpy array, with dimensions
        [bid, bax, datatype (time or value), timepoints]

    Returns
    -------
    tuple of two np.arrays
        The first array contains the fitted k values, the second the fitted
        fmax values. Both are of dimension [bid concs, bax_concs].
    """

    fmax_data = np.zeros((len(bid_concs), len(bax_concs)))
    k_data = np.zeros((len(bid_concs), len(bax_concs)))

    for bid_ix in range(data_norm.shape[0]):
        for bax_ix in range(data_norm.shape[1]):

            t = data_norm[bid_ix, bax_ix, TIME, :]
            v = data_norm[bid_ix, bax_ix, VALUE, :]

            # Get an instance of the fitting function
            fit = tf.OneExpFmax()
            # Subtract 1 so that the curves start from 0
            (k, fmax) = fit.fit_timecourse(t, v - 1)
            k_data[bid_ix, bax_ix] = k
            fmax_data[bid_ix, bax_ix] = fmax

            # Some diagnostic plots (commented out)
            #plt.figure()
            #plt.plot(t, v, 'k')
            #plt.plot(t[:numpts], intercept + lin_fit[0] * t[:numpts], 'r')

            #plt.figure()
            #plt.plot(t, v / intercept - 1)
            #plt.plot(t, fit.fit_func(t, [k, fmax]))

    return (k_data, fmax_data)

def plot_bax_titration_timecourses(data_matrix, bid_ix, k_data, fmax_data,
                                   plot_filename=None):
    """Plot F/F0 vs. time for a Bid conc and a range of Bax concs."""
    fig = plt.figure(figsize=(1.5, 1.5), dpi=300)
    plt.subplots_adjust(left=0.24, bottom=0.21)
    ax = fig.gca()
    for bax_ix in range(data_matrix.shape[1]):
        if bax_ix not in [0, 6, 10]:
            continue
        t = data_matrix[bid_ix, bax_ix, TIME, :]
        v = data_matrix[bid_ix, bax_ix, VALUE, :]
        ax.plot(t, v, alpha=0.9, linewidth=0.5)
        # Plot exponential fits
        k = k_data[bid_ix, bax_ix]
        fmax = fmax_data[bid_ix, bax_ix]
        fit = tf.OneExpFmax()
        ax.plot(t, fit.fit_func(t, [k, fmax]) + 1, 'k')

    ax.set_xticks(np.linspace(0, 2.5e4, 6))
    ax.set_xticklabels([int(f) for f in np.linspace(0, 25, 6)])
    ax.set_ylabel('$F/F_0$')
    ax.set_xlabel(r'Time (sec $\times 10^3$)')
    format_axis(ax)

    if plot_filename:
        fig.savefig('%s.pdf' % plot_filename)
        fig.savefig('%s.png' % plot_filename, dpi=300)

def plot_k_fmax_scaling(k_data, fmax_data, bid_ix, bax_concs,
                        plot_filename=None):
    """Plot Fmax and k vs. Bax concentration for a given Bid concentration."""

    assert k_data.shape == fmax_data.shape, \
           "k_data and f_max data must have same dimensions"

    fig = plt.figure('exp_fits_k_fmax_var', figsize=(1.9, 1.5), dpi=300)
    ax1 = fig.gca()
    # Plot Fmax vs. concentration on the left-hand axis
    ax1.plot(bax_concs, fmax_data[bid_ix, :], marker='o', markersize=3,
             color='b')
    ax1.set_xlabel('[Bax] (nM)')
    ax1.set_ylabel('$F_{max}$', color='b')
    ax1.set_xlim([0, 1100])
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
    #ax1.set_yscale('log')
    # Plot k vs. concentration on the right-hand axis
    ax2 = ax1.twinx()
    #ax2.set_xlim([0.05, 30])
    ax2.plot(bax_concs, k_data[bid_ix,:], marker='o', markersize=3,
             color='r')
    ax2.set_ylabel(r'k (sec$^{-1} \times\ 10^{-5}$)', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    ax2.set_xlim([0, 1100])
    ax2.set_ylim(5e-5, 3e-4)
    ax2.set_yticks(np.linspace(5e-5, 3e-4, 6))
    ax2.set_yticklabels([str(int(f)) for f in np.linspace(5, 30, 6)])
    format_axis(ax1)
    format_axis(ax2, yticks_position='right')
    plt.subplots_adjust(left=0.18, bottom=0.19, right=0.75)

    if plot_filename:
        fig.savefig('%s.pdf' % plot_filename)
        fig.savefig('%s.png' % plot_filename, dpi=300)

def plot_endpoints_vs_bax(data_norm, time_pts, bid_ix, bax_concs,
                          plot_filename=None, avg_pts=10):
    # Plot endpoints at 2, 3, and 5 hours
    fig = plt.figure(figsize=(1.5, 1.5), dpi=300)
    ax = fig.gca()

    for time_ix in time_pts:
        # Get endpoint vector
        start_ix = time_ix - avg_pts
        f_avg = np.mean(data_norm[bid_ix, :, VALUE, start_ix:time_ix], axis=1)
        f_sd = np.std(data_norm[bid_ix, :, VALUE, start_ix:time_ix], axis=1,
                      ddof=1)
        ax.errorbar(bax_concs, f_avg, yerr=f_sd)

    ax.set_ylabel('NBD $F/F_0$')
    ax.set_xlabel('[Bax] (nM)')
    format_axis(ax)
    plt.subplots_adjust(left=0.24, bottom=0.21)
    ax.set_xlim([10, 1500])
    ax.set_xscale('log')

    if plot_filename:
        fig.savefig('%s.pdf' % plot_filename)
        fig.savefig('%s.png' % plot_filename, dpi=300)

if __name__ == '__main__':
    #plot_data()
    set_fig_params_for_publication()

    # Fit the data with exponential nctions
    (k_data, fmax_data) = calc_exp_fits(data_norm)

    # Plot raw timecourses with fits
    plot_bax_titration_timecourses(data_norm, 4, k_data, fmax_data,
                plot_filename='141119_Bid_20nm_timecourses')

    # Plot k and fmax scaling
    plot_k_fmax_scaling(k_data, fmax_data, 4, bax_concs,
                               plot_filename='141119_Bid_20nm_scaling')

    # Plot the F/F0 values at various pts vs. Bax concentration
    plot_endpoints_vs_bax(data_norm, [125, 250, 500], 4, bax_concs,
                          plot_filename='141119_Bid_20nm_endpts',
                          avg_pts=10)

    sys.exit()


    plt.figure()
    plt.errorbar(bax_concs, endpt_data[4, :], yerr=endpt_sd_data[4, :],
                 marker='o')


    set_fig_params_for_publication()

    k_data = []

    plt.figure('k', figsize=(1.5, 1.5), dpi=300)
    plt.subplots_adjust(left=0.24, bottom=0.21)
    plt.xlabel('[Total Bax] (nM)')
    plt.ylabel(r'k (sec$^{-1} \times 10^{-5}$)')
    ax = plt.gca()
    ax.set_ylim([5e-5, 37e-5])
    ax.set_xlim([10, 2000])
    ax.set_yticks(np.linspace(5e-5, 35e-5, 7))
    ax.set_yticklabels([int(f) for f in np.linspace(5, 35, 7)])
    ax.set_xscale('log')
    format_axis(ax)

    plt.figure('fmax', figsize=(1.5, 1.5), dpi=300)
    plt.subplots_adjust(left=0.24, bottom=0.21)
    plt.xlabel('[Total Bax] (nM)')
    plt.ylabel(r'$F_{max}$')
    ax = plt.gca()
    ax.set_ylim([1, 5])
    ax.set_xlim([10, 2000])
    #ax.set_yticks(np.linspace(5e-5, 35e-5, 7))
    #ax.set_yticklabels([int(f) for f in np.linspace(5, 35, 7)])
    ax.set_xscale('log')
    format_axis(ax)


import sys
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import linregress

from tbidbaxlipo.util.plate_assay import plot_all, TIME, VALUE
from tbidbaxlipo.util import fitting, colors, set_fig_params_for_publication, \
                             format_axis, fontsize
from preprocess_data import bim_bh3, bid_80, bid_40, bid_20, bid_10, \
                            bid_5, bid_2, bid_0, bax_concs, bid_concs, \
                            bg_averages, bg_diff, data_matrix, data_norm
from tbidbaxlipo.plots import titration_fits as tf

#plt.ion()

brew_colors11 = [
    '#a50026',
    '#d73027',
    '#f46d43',
    '#fdae61',
    '#fee090',
    '#ffffbf',
    '#e0f3f8',
    '#abd9e9',
    '#74add1',
    '#4575b4',
    '#313695',
]
brew_colors11.reverse()

brew_colors7 = [
    '#d73027',
    '#fc8d59',
    '#fee090',
    '#ffffbf',
    '#e0f3f8',
    '#91bfdb',
    '#4575b4',
]
brew_colors7.reverse()

brew_colors6 = [
    '#d73027',
    '#fc8d59',
    '#fee090',
    '#e0f3f8',
    '#91bfdb',
    '#4575b4',
]
brew_colors6.reverse()

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


def plot_mm(k_data, bid_concs, bax_concs, plot_filename=None):
    # Create the figure
    fig = plt.figure('mm fit', figsize=(1.5, 1.5), dpi=300)
    ax = fig.gca()

    km_list = []
    kcat_list = []
    v0_list = []

    for bid_ix in range(1, len(bid_concs)):
        bid_conc = bid_concs[bid_ix]
        # Fitting parameters
        kcat = fitting.Parameter(0.06)
        km = fitting.Parameter(250.)
        v0 = fitting.Parameter(5e-5)

        # Fit func for this bid concentration
        def fit_func(bax_concs):
            return ((kcat() * bid_conc) / (km() + bax_concs)) + v0()

        # Titration k data for this Bid concentration
        bid_k = k_data[bid_ix, :]
        # Fit to the MM equation
        fitting.fit(fit_func, [kcat, km, v0], bid_k, np.array(bax_concs))

        # Save the fitted values
        km_list.append(km())
        kcat_list.append(kcat())
        v0_list.append(v0())

        # Plotting
        #c = colors[bid_ix]
        c = brew_colors6[bid_ix-1]
        bid_str = '2.5' if bid_conc == 2.5 else str(int(bid_conc))
        # Plot the MM fit
        ax.plot(bax_concs, fit_func(bax_concs), color=c,
                 label='%s nM' % bid_str)
        # Plot the data
        ax.plot(bax_concs, bid_k, marker='o', color=c,
                linestyle='', markersize=3)

    # Format the plot
    plt.subplots_adjust(left=0.21, bottom=0.19)
    ax.set_xlabel('[Total Bax] (nM)')
    ax.set_ylabel(r'k (sec$^{-1} \times 10^{-5}$)')
    ax.set_ylim([-1e-5, 37e-5])
    ax.set_xlim([10, 2000])
    ax.set_yticks(np.linspace(0, 35e-5, 8))
    ax.set_yticklabels([int(f) for f in np.linspace(0, 35, 8)])
    ax.set_xscale('log')
    format_axis(ax)

    #max_rate_list = []
    #for bid_ix, bid_conc in enumerate(bid_concs):
    #    km = km_list[bid_ix]
    #    kcat = kcat_list[bid_ix]
    #    v0 = v0_list[bid_ix]
    #   #max_rate = ((kcat * bid_conc) / km) + v0
    #    max_rate = (kcat * bid_conc) + v0
    #    max_rate_list.append(max_rate)

    kcat_fig = plt.figure('kcat', figsize=(1.5, 1.5), dpi=300)
    kcat_ax = kcat_fig.gca()
    #kcat_ax.plot(bid_concs[1:], max_rate_list[1:], marker='o')
    max_rate_list = [k_data[i, 0] for i in range(len(bid_concs))]
    kcat_ax.plot(bid_concs, max_rate_list, marker='o', markersize=3,
                 label='Observed $k$')

    #kcat_ax.set_xscale('log')
    kcat_ax.set_xlabel('[cBid] (nM)')
    kcat_ax.set_ylabel('Rate')
    kcat_ax.set_ylabel(r'k (sec$^{-1} \times 10^{-5}$)')
    kcat_ax.set_ylim([-1e-5, 37e-5])
    kcat_ax.set_yticks(np.linspace(0, 35e-5, 8))
    kcat_ax.set_yticklabels([int(f) for f in np.linspace(0, 35, 8)])
    kcat_ax.set_xlim(-2, 82)
    # Define line connecting 0 Bid and 2.5 Bid concs
    bid_slope = (max_rate_list[1] - max_rate_list[0]) / float(bid_concs[1])
    kcat_ax.plot(bid_concs, bid_slope * bid_concs + max_rate_list[0],
                 color='r', label='Fixed $k_{cat}$')
    plt.legend(loc='lower right', frameon=False, fontsize=fontsize)
    kcat_fig.subplots_adjust(left=0.22, bottom=0.2)
    format_axis(kcat_ax)

    if plot_filename:
        fig.savefig('%s.pdf' % plot_filename)
        fig.savefig('%s.png' % plot_filename, dpi=300)
        kcat_fig.savefig('%s_kcat.pdf' % plot_filename)
        kcat_fig.savefig('%s_kcat.png' % plot_filename, dpi=300)

    #plt.plot(bid_concs[1:], km_list[1:], marker='o')
    #plt.figure('Kcat')
    #plt.plot(bid_concs[1:], kcat_list[1:], marker='o')
    #plt.figure('v0')
    #plt.plot(bid_concs[1:], v0_list[1:], marker='o')

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
    fmax_sd_data = np.zeros((len(bid_concs), len(bax_concs)))
    k_sd_data = np.zeros((len(bid_concs), len(bax_concs)))

    for bid_ix in range(data_norm.shape[0]):
        for bax_ix in range(data_norm.shape[1]):

            t = data_norm[bid_ix, bax_ix, TIME, :]
            v = data_norm[bid_ix, bax_ix, VALUE, :]

            # Get an instance of the fitting function
            k = fitting.Parameter(1e-4)
            fmax = fitting.Parameter(3.)

            def fit_func(t):
                return 1 + fmax() * (1 - np.exp(-k()*t))
            # Subtract 1 so that the curves start from 0
            (residuals, result) = fitting.fit(fit_func, [k, fmax], v, t)
            cov_matrix = result[1]
            #(k, fmax) = fit.fit_timecourse(t, v)
            k_data[bid_ix, bax_ix] = k()
            fmax_data[bid_ix, bax_ix] = fmax()

            #import ipdb; ipdb.set_trace()
            # Get variance estimates
            [k_sd, fmax_sd] = np.sqrt(np.diag(np.var(residuals, ddof=1) *
                                              cov_matrix))
            k_sd_data[bid_ix, bax_ix] = k_sd
            fmax_sd_data[bid_ix, bax_ix] = fmax_sd

            # Some diagnostic plots (commented out)
            #plt.figure()
            #plt.plot(t, v, 'k')
            #plt.plot(t[:numpts], intercept + lin_fit[0] * t[:numpts], 'r')

            #plt.figure()
            #plt.plot(t, v / intercept - 1)
            #plt.plot(t, fit.fit_func(t, [k, fmax]))

    return (k_data, k_sd_data, fmax_data, fmax_sd_data)

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
        ax.plot(t, v, alpha=0.9, linewidth=1,
                label='%d nM Bax' % bax_concs[bax_ix])
        # Plot exponential fits
        k = k_data[bid_ix, bax_ix]
        fmax = fmax_data[bid_ix, bax_ix]
        fit = tf.OneExpFmax()
        ax.plot(t, fit.fit_func(t, [k, fmax]) + 1, 'k')
    plt.legend(loc='lower right', fontsize=fontsize, frameon=False)
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
    ax1.set_xscale('log')
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
    colors = ['r', 'g', 'b']
    for i, time_ix in enumerate(time_pts):
        # Get endpoint vector
        start_ix = time_ix - avg_pts
        f_avg = np.mean(data_norm[bid_ix, :, VALUE, start_ix:time_ix], axis=1)
        f_sd = np.std(data_norm[bid_ix, :, VALUE, start_ix:time_ix], axis=1,
                      ddof=1)
        ax.errorbar(bax_concs, f_avg, yerr=f_sd, color=colors[i])


    ax.set_ylabel('NBD $F/F_0$')
    ax.set_xlabel('[Bax] (nM)')
    format_axis(ax)
    plt.subplots_adjust(left=0.24, bottom=0.21)
    ax.set_xlim([10, 1500])
    ax.set_xscale('log')

    if plot_filename:
        fig.savefig('%s.pdf' % plot_filename)
        fig.savefig('%s.png' % plot_filename, dpi=300)

def plot_k_data(k_data, k_sd_data, bid_concs, bax_concs, plot_filename=None):
    fig = plt.figure(figsize=(1.5, 1.5), dpi=300)
    ax = fig.gca()

    for bid_ix in reversed(range(len(bid_concs))):
        bid_conc = bid_concs[bid_ix]
        # Titration k data for this Bid concentration
        bid_k = k_data[bid_ix, :]
        bid_k_sd = k_sd_data[bid_ix, :]
        bid_k_ub = 10 ** (np.log10(bid_k) + bid_k_sd) - bid_k
        bid_k_lb = bid_k - 10 ** (np.log10(bid_k) - bid_k_sd)
        #import ipdb; ipdb.set_trace()
        # Plotting
        c = brew_colors6[bid_ix]
        bid_str = '2.5' if bid_conc == 2.5 else str(int(bid_conc))
        # Plot the data
        ax.errorbar(bax_concs, bid_k, yerr=[bid_k_lb, bid_k_ub], marker='o',
                    color=c, linestyle='-', markersize=3)

    # Format the plot
    plt.subplots_adjust(left=0.21, bottom=0.19)
    ax.set_xlabel('[Total Bax] (nM)')
    ax.set_ylabel(r'$k$ (sec$^{-1} \times 10^{-5}$)')
    ax.set_ylim([-1e-5, 37e-5])
    ax.set_xlim([10, 2000])
    ax.set_yticks(np.linspace(0, 35e-5, 8))
    ax.set_yticklabels([int(f) for f in np.linspace(0, 35, 8)])
    ax.set_xscale('log')
    format_axis(ax)

    if plot_filename:
        fig.savefig('%s.pdf' % plot_filename)
        fig.savefig('%s.png' % plot_filename, dpi=300)

    # Now, plot the same data looking at the scaling with Bid
    bid_fig = plt.figure(figsize=(1.5, 1.5), dpi=300)
    ax = bid_fig.gca()

    for bax_ix in reversed(range(len(bax_concs))):
        bax_conc = bax_concs[bax_ix]
        # Titration k data for this Bax concentration
        bax_k = k_data[:, bax_ix]
        bax_k_sd = k_sd_data[:, bax_ix]
        bax_k_ub = 10 ** (np.log10(bax_k) + bax_k_sd) - bax_k
        bax_k_lb = bax_k - 10 ** (np.log10(bax_k) - bax_k_sd)
        #import ipdb; ipdb.set_trace()
        # Plotting
        c = brew_colors11[bax_ix]
        # Plot the data
        ax.errorbar(bid_concs, bax_k, yerr=[bax_k_lb, bax_k_ub], marker='o',
                    color=c, linestyle='-', markersize=3, capsize=3)

    # Format the plot
    plt.subplots_adjust(left=0.21, bottom=0.19)
    ax.set_xlabel('[cBid] (nM)')
    ax.set_ylabel(r'$k$ (sec$^{-1} \times 10^{-5}$)')
    ax.set_ylim([-1e-5, 37e-5])
    ax.set_xlim([-2, 82])
    ax.set_yticks(np.linspace(0, 35e-5, 8))
    ax.set_yticklabels([int(f) for f in np.linspace(0, 35, 8)])
    format_axis(ax)

    if plot_filename:
        bid_fig.savefig('%s_bid.pdf' % plot_filename)
        bid_fig.savefig('%s_bid.png' % plot_filename, dpi=300)

def plot_fmax_data(fmax_data, bid_concs, bax_concs, plot_filename=None):
    fig = plt.figure(figsize=(1.5, 1.5), dpi=300)
    ax = fig.gca()

    for bid_ix, bid_conc in enumerate(bid_concs):
        # Titration fmax data for this Bid concentration
        bid_fmax = fmax_data[bid_ix, :]

        # Plotting
        #c = colors[bid_ix]
        c = brew_colors6[bid_ix]
        bid_str = '2.5' if bid_conc == 2.5 else str(int(bid_conc))
        # Plot the data
        ax.plot(bax_concs, bid_fmax, marker='o', color=c,
                 linestyle='-', markersize=3)

    # Format the plot
    plt.subplots_adjust(left=0.23, bottom=0.19, right=0.92)
    ax.set_xlabel('[Total Bax] (nM)')
    ax.set_ylabel(r'$F_{max}$ (NBD F/F$_0$)')
    #ax.set_ylim([-1e-5, 37e-5])
    ax.set_xlim([10, 2000])
    # ax.set_yticks(np.linspace(0, 35e-5, 8))
    # ax.set_yticklabels([int(f) for f in np.linspace(0, 35, 8)])
    ax.set_xscale('log')
    format_axis(ax)

    if plot_filename:
        fig.savefig('%s.pdf' % plot_filename)
        fig.savefig('%s.png' % plot_filename, dpi=300)

if __name__ == '__main__':
    #plot_data()
    set_fig_params_for_publication()
    #plt.ion()
    # Fit the data with exponential nctions
    (k_data, k_sd_data, fmax_data, fmax_sd_data) = calc_exp_fits(data_norm)

    # Plot the scaling of k vs. [Bax] for all [Bid]
    plot_k_data(k_data[1:], k_sd_data[1:], bid_concs[1:], bax_concs,
                plot_filename='141119_k_scaling')

    # Plot the scaling of k vs. [Bax] for all [Bid]
    #plot_k_data(k_data[0:], k_sd_data[0:], bid_concs[0:], bax_concs,
    #            plot_filename='141119_k_scaling')
    #sys.exit()

    # Plot the scaling of fmax vs. [Bax] for all [Bid]
    plot_fmax_data(fmax_data[1:], bid_concs[1:], bax_concs,
                   plot_filename='141119_fmax_scaling')

    # Plot k scaling with Michaelis-Menten local fits
    plot_mm(k_data, bid_concs, bax_concs,
            plot_filename='141119_k_mm_fits')

    # Plot raw timecourses with fits
    plot_bax_titration_timecourses(data_norm, 1, k_data, fmax_data,
                plot_filename='141119_Bid_2nm_timecourses')

    # Plot raw timecourses with fits
    plot_bax_titration_timecourses(data_norm, 4, k_data, fmax_data,
                plot_filename='141119_Bid_20nm_timecourses')


    # Plot raw timecourses with fits
    plot_bax_titration_timecourses(data_norm, 6, k_data, fmax_data,
                plot_filename='141119_Bid_80nm_timecourses')

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


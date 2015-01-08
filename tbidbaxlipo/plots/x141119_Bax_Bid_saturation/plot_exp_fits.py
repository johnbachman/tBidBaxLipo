import sys
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import linregress

from tbidbaxlipo.util.plate_assay import plot_all, TIME, VALUE
from tbidbaxlipo.util import fitting, colors, set_fig_params_for_publication, \
                             format_axis
from preprocess_data import bim_bh3, bid_80, bid_40, bid_20, bid_10, \
                            bid_5, bid_2, bid_0, bax_concs
from tbidbaxlipo.plots import titration_fits as tf

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

if __name__ == '__main__':
    plt.ion()
    #plot_data()

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

    for bid in [bid_2, bid_5, bid_10, bid_20, bid_40, bid_80]:
        fmax_list = []
        k_list = []

        for well_name in bid.keys():
            well = bid[well_name]
            t = well[TIME]
            v = well[VALUE]
            #figure()
            # Do linear regression on the first 20 points, normalize to intcpt
            numpts = 20
            lin_fit = linregress(t[:numpts], v[:numpts])
            intercept = lin_fit[1]
            #plot(t, v, 'k')
            #plot(t[:numpts], intercept + lin_fit[0] * t[:numpts], 'r')
            fit = tf.OneExpFmax()
            (k, fmax) = fit.fit_timecourse(t, v / intercept - 1)
            k_list.append(k)
            fmax_list.append(fmax)
            #figure()
            #plot(t, v / intercept - 1)
            #plot(t, fit.fit_func(t, [k, fmax]))

        k_data.append(k_list)

        #figure('k')
        #plot(bax_concs, k_list, marker='o')
        #title('k')
        plt.figure('k')
        plt.plot(bax_concs, k_list, marker='o', markersize=3)

        plt.figure('fmax')
        plt.plot(bax_concs, fmax_list, marker='o', markersize=3)
    ax = plt.gca()
    format_axis(ax)

    bid_concs = np.array([2.5, 5., 10., 20., 40., 80.])
    plot_mm(k_data, bax_concs, bid_concs)


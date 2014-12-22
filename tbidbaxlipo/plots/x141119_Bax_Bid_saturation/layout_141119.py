from tbidbaxlipo.util.plate_assay import *
import itertools
import pickle
import collections
import sys
import os
import tbidbaxlipo.data
from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import fitting, set_fig_params_for_publication, emcee_fit
from tbidbaxlipo.plots.titration_fits import TwoExp, OneExpFmax
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from tbidbaxlipo.models.nbd import multiconf
from tbidbaxlipo.models import one_cpt
from scipy.stats import linregress
from tbidbaxlipo.plots import titration_fits as tf
from tbidbaxlipo.util import fitting, colors, set_fig_params_for_publication, \
                             format_axis

# All wells with liposomes ~0.15 mg/mL

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                  '141119_Bax_Bid_saturation_timecourse.txt'))

bax_concs = np.array([1000., 500., 250., 125., 62.5, 31.2, 15.6, 7.8,
                      3.9, 2.0, 0.]) + 25.

# Timecourse wells
timecourse_wells = read_flexstation_kinetics(timecourse_file)
"""The raw (unnormalized) timecourses."""

def time_offset_vector():
    read_time = 320 # seconds (5:20)
    row_times = {'A': 0,
                 'B': 20,
                 'C': 37,
                 'D': 60,
                 'E': 80,
                 'F': 100,
                 'G': 130,
                 'H': 150,
             }
    offset_dict = {}
    # New time vector is TIME + offset, where
    # offset = (read_time - row_time)
    for row in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
        for col in range(1, 13):
            well = '%s%s' % (row, col)
            offset_dict[well] = read_time - row_times[row]
    return offset_dict

def add_offset_vector(timecourse_wells, offset_vector):
    for well_name in timecourse_wells.keys():
        well = timecourse_wells[well_name]
        well[TIME] += offset_vector[well_name]

def row_wells(row_name, max_col=12):
    return ['%s%s' % (row_name, i) for i in range(1, max_col+1)]

# Should update time vector of timecourse wells in place
add_offset_vector(timecourse_wells, time_offset_vector())

no_nbd_lipo_well_names = ['B12', 'D12', 'F12'] # Ignoring H12 as an outlier
no_nbd_lipo_wells = extract(no_nbd_lipo_well_names, timecourse_wells)

no_lipo_well_names = ['A12', 'C12', 'E12', 'G12']
no_lipo_wells = extract(no_lipo_well_names, timecourse_wells)

bg_layout = collections.OrderedDict([
         ('No lipos', ['A12', 'C12', 'E12', 'G12']),
         ('No NBD or lipos', ['B12', 'D12', 'F12']),
         ('NBD and lipos', ['G9', 'G10', 'G11']),
        ])
(bg_averages, bg_std) = averages(timecourse_wells, bg_layout)

# Difference between NBD-Bax and background (buffer only) gives fluorescence
# of NBD-Bax alone
bg_diff = {}
bg_diff['BG diff'] = []
bg_diff['BG diff'].append(bg_averages['No lipos'][TIME])
bg_diff['BG diff'].append(bg_averages['No lipos'][VALUE] -
                          bg_averages['No NBD or lipos'][VALUE])
bg_diff['BG diff2'] = []
bg_diff['BG diff2'].append(bg_averages['NBD and lipos'][TIME])
bg_diff['BG diff2'].append(bg_averages['NBD and lipos'][VALUE] -
                           bg_averages['No lipos'][VALUE])
bg_diff['BG diff3'] = []
bg_diff['BG diff3'].append(bg_averages['No NBD or lipos'][TIME])
bg_diff['BG diff3'].append(bg_averages['No NBD or lipos'][VALUE] + 2.5)

# So liposomes alone has fluorescence of roughly 2.55
# So need to subtract background and 2.55 from each curve to get
# baseline fluorescence. The other curves are normalized against this.

bg_to_subtract = bg_averages['No NBD or lipos'][VALUE] + 2.5

bg_sub_wells = subtract_background(timecourse_wells, bg_to_subtract)

bid_80 = extract(row_wells('A', 11), bg_sub_wells)
bid_40 = extract(row_wells('B', 11), bg_sub_wells)
bid_20 = extract(row_wells('C', 11), bg_sub_wells)
bid_10 = extract(row_wells('D', 11), bg_sub_wells)
bid_5 = extract(row_wells('E', 11), bg_sub_wells)
bid_2 = extract(row_wells('F', 11), bg_sub_wells)
bid_0 = extract(row_wells('G', 11), bg_sub_wells)
bim_bh3 = extract(row_wells('H', 11), bg_sub_wells)

def plot_data():
    """Plots the data and various transformations of it."""
    """
    figure()
    plot_all(no_nbd_lipo_wells)
    title("No NBD-Bax or Liposomes")

    # Lipo bg wells
    figure()
    plot_all(no_lipo_wells)
    title("No liposomes")

    # Lipo bg wells
    figure()
    plot_all(bg_averages, errors=bg_std)
    title("Background")

    # BG difference
    figure()
    plot_all(bg_diff)

    figure()
    plot_all(bg_sub_wells)
    """

    figure()
    plot_all(bim_bh3)
    title('Bim BH3')

    figure()
    plot_all(bid_80)
    title('Bid 80nM')

    figure()
    plot_all(bid_40)
    title('Bid 40 nM')

    figure()
    plot_all(bid_20)
    title('Bid 20 nM')

    figure()
    plot_all(bid_10)
    title('Bid 10 nM')

    figure()
    plot_all(bid_5)
    title('Bid 5 nM')

    figure()
    plot_all(bid_2)
    title('Bid 2 nM')

    figure()
    plot_all(bid_0)
    title('Bid 0 nM')

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


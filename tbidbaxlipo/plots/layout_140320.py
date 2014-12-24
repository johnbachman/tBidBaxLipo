from tbidbaxlipo.util.plate_assay import *
import itertools
import pickle
import collections
import sys
import os
import tbidbaxlipo.data
from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import fitting, format_axis, \
                             set_fig_params_for_publication
from pysb import *

layout = collections.OrderedDict([
        ('Bax 590 nM, NBD-Bax 96 nM',  ['A1', 'B1']),
        ('Bax 393 nM, NBD-Bax 96 nM',  ['A2', 'B2']),
        ('Bax 262 nM, NBD-Bax 96 nM',  ['A3', 'B3']),
        ('Bax 175 nM, NBD-Bax 96 nM',  ['A4', 'B4']),
        ('Bax 116 nM, NBD-Bax 96 nM',  ['A5', 'B5']),
        ('Bax 78 nM, NBD-Bax 96 nM',  ['A6', 'B6']),
        ('Bax 52 nM, NBD-Bax 96 nM',  ['A7', 'B7']),
        ('Bax 35 nM, NBD-Bax 96 nM',  ['A8', 'B8']),
        ('Bax 23 nM, NBD-Bax 96 nM',  ['A9', 'B9']),
        ('Bax 15 nM, NBD-Bax 96 nM',  ['A10', 'B10']),
        ('Bax 10 nM, NBD-Bax 96 nM',  ['A11', 'B11']),
        ('Bax 0 nM, NBD-Bax 96 nM',  ['A12', 'B12']),
        ('Bax 590 nM',  ['C1']),
        ('Bax 393 nM',  ['C2']),
        ('Bax 262 nM',  ['C3']),
        ('Bax 175 nM',  ['C4']),
        ('Bax 116 nM',  ['C5']),
        ('Bax 78 nM',  ['C6']),
        ('Bax 52 nM',  ['C7']),
        ('Bax 35 nM',  ['C8']),
        ('Bax 23 nM',  ['C9']),
        ('Bax 15 nM',  ['C10']),
        ('Bax 10 nM',  ['C11']),
        ('Bax 0 nM',  ['C12']),
    ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '140320_NBD_Bax_BimBH3_unlab_Bax_titration.txt'))

# Assemble all the wells included in the layout
# http://stackoverflow.com/questions/406121/
# flattening-a-shallow-list-in-python
wells_to_read = list(itertools.chain(*layout.values()))

# Timecourse wells
timecourse_wells = read_flexstation_kinetics(timecourse_file)
"""The raw (unnormalized) timecourses."""
timecourse_wells = extract(wells_to_read, timecourse_wells)

# Averaged
(timecourse_averages, timecourse_stds) = averages(timecourse_wells, layout)

bax_bg_conditions = [
        'Bax 590 nM',
        'Bax 393 nM',
        'Bax 262 nM',
        'Bax 175 nM',
        'Bax 116 nM',
        'Bax 78 nM',
        'Bax 52 nM',
        'Bax 35 nM',
        'Bax 23 nM',
        'Bax 15 nM',
        'Bax 10 nM',
        'Bax 0 nM']

bax_bg_layout = extract(bax_bg_conditions, layout)
bax_bg_wells = extract([layout[cond][0] for cond in bax_bg_conditions],
                        timecourse_wells)

nbd_conditions = [
        'Bax 590 nM, NBD-Bax 96 nM',
        'Bax 393 nM, NBD-Bax 96 nM',
        'Bax 262 nM, NBD-Bax 96 nM',
        'Bax 175 nM, NBD-Bax 96 nM',
        'Bax 116 nM, NBD-Bax 96 nM',
        'Bax 78 nM, NBD-Bax 96 nM',
        'Bax 52 nM, NBD-Bax 96 nM',
        'Bax 35 nM, NBD-Bax 96 nM',
        'Bax 23 nM, NBD-Bax 96 nM',
        'Bax 15 nM, NBD-Bax 96 nM',
        'Bax 10 nM, NBD-Bax 96 nM',
        'Bax 0 nM, NBD-Bax 96 nM',]

nbd_layout = extract(nbd_conditions, layout)
nbd_wells = extract(nbd_conditions, timecourse_averages)
#nbd_wells = extract([layout[cond][0] for cond in nbd_conditions],
#                         timecourse_wells)

# Background subtracted
bgsub_wells = subtract_background_set(nbd_wells, bax_bg_wells)

# Background subtracted, averaged
#(bgsub_averages, bgsub_stds) = averages(bgsub_wells, layout)

# First timepoint shifted to 0 (better for fitting)
#reset_bgsub_means = reset_first_timepoint_to_zero(bgsub_norm_averages)
#Timecourses normalized, BG-subtracted, averaged, then with first point
#shifted to t = 0.
#reset_bgsub_sds = reset_first_timepoint_to_zero(bgsub_norm_stds)

def plot_data():
    """Plots the data and various transformations of it."""
    ion()

    # Timecourse wells
    figure()
    plot_all(timecourse_wells)
    title("Raw timecourses")

    # NBD wells
    figure()
    plot_all(nbd_wells)
    title("NBD-Bax wells, Raw")

    # Lipo bg wells
    figure()
    plot_all(bax_bg_wells)
    title("Unlabeled Bax background wells, Raw")

    # Averages
    figure()
    plot_all(timecourse_averages)
    title("Average timecourses, all")

    # Background-subtracted wells
    figure()
    plot_all(bgsub_wells)
    title("Background-subtracted wells")

    # Background-subtracted wells, averaged
    #figure()
    #plot_all(bgsub_averages)
    #title("Background-subtracted wells")

def plot_exp_fits():
    fmax_list = []
    k1_list = []
    conc_list = []

    for conc_name in bgsub_wells.keys():
        conc_list.append(float(conc_name.split(' ')[1]))
        # Try fitting the high conc trajectory
        y = bgsub_wells[conc_name][VALUE]
        y = y / np.mean(y[0:2])
        t = bgsub_wells[conc_name][TIME]
        fmax = fitting.Parameter(25.)
        k1 = fitting.Parameter(np.log(2)/2000.)
        k2 = fitting.Parameter(1e-3)
        b = fitting.Parameter(np.mean(y[0:2]))
        k_bleach = fitting.Parameter(1.17e-5)
        #def exp_func(t):
        #    return (b() + fmax()*(1 - np.exp(-k1() *
        #            (1 - np.exp(-k2()*t))*t))) * np.exp(-k_bleach()*t)
        # One-parameter exp
        def exp_func(t):
            return (b() + fmax()*(1 - np.exp(-k1()*t))) * \
                    np.exp(-k_bleach()*t)

        #fitting.fit(exp_func, [fmax, k1, k2], y, t)
        fitting.fit(exp_func, [fmax, k1], y, t)

        # Add to list
        fmax_list.append(fmax())
        k1_list.append(k1())

        #plt.figure()
        #plt.plot(t, y, marker='o', linestyle='')
        #plt.plot(t, exp_func(t), color='r', linewidth='2')
        #plt.title(conc_name)

    set_fig_params_for_publication()

    plt.figure('exp_fits', figsize=(1.7, 1.5), dpi=300)
    plt.plot(conc_list[:-1], fmax_list[:-1], marker='o', markersize=3,
             color='b')
    plt.xlabel('[Bax] (nM)')
    ax1 = plt.gca()
    ax1.set_xscale('log')
    ax1.set_xlim([8, 1000])
    ax1.set_ylabel('$F_{max}$', color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')

    ax2 = ax1.twinx()
    ax2.set_xlim([5, 1000])
    ax2.plot(conc_list[:-1], k1_list[:-1], marker='o', markersize=3, color='r')
    ax2.set_ylabel(r'k (sec $\times$ 10^{-3})', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    ax2.set_yticks(np.linspace(6.6e-4, 7.8e-4, 7))
    ax2.set_yticklabels(['%.1f' % f for f in np.linspace(6.6, 7.8, 7)])

    format_axis(ax1)
    format_axis(ax2, yticks_position='right')
    plt.subplots_adjust(left=0.20, bottom=0.19, right=0.80)

if __name__ == '__main__':
    plt.ion()
    # plot_data()

    plot_exp_fits()
    plt.savefig('140320_exp_fits.pdf')

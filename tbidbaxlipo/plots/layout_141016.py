from tbidbaxlipo.util.plate_assay import *
import itertools
import pickle
import collections
import sys
import os
import tbidbaxlipo.data
from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.plots import titration_fits
import pymc
from scipy import stats
from tbidbaxlipo.util import fig_orange, fig_purple, \
                             set_fig_params_for_publication, fontsize, format_axis
from copy import copy
import matplotlib

# All wells have 100 nM cBid and 100 nM Bax
layout = collections.OrderedDict([
        ('Lipos 21.7 nM, Pre-inc',  ['C1', 'D1']),
        ('Lipos 10.9 nM, Pre-inc',  ['C2', 'D2']),
        ('Lipos 5.4 nM, Pre-inc',  ['C3', 'D3']),
        ('Lipos 2.7 nM, Pre-inc',  ['C4', 'D4']),
        ('Lipos 1.4 nM, Pre-inc',  ['C5', 'D5']),
        ('Lipos 0.68 nM, Pre-inc',  ['C6', 'D6']),
        ('Lipos 0.34 nM, Pre-inc',  ['C7', 'D7']),
        ('Lipos 0.17 nM, Pre-inc',  ['C8', 'D8']),
        ('Lipos 0.085 nM, Pre-inc',  ['C9', 'D9']),
        ('Lipos 0.042 nM, Pre-inc',  ['C10', 'D10']),
        ('Lipos 0.021 nM, Pre-inc',  ['C11', 'D11']),
        ('Lipos 0 nM, Pre-inc',  ['C12', 'D12']),
        ('Lipos 21.7 nM',  ['E1', 'F1']),
        ('Lipos 10.9 nM',  ['E2', 'F2']),
        ('Lipos 5.4 nM',  ['E3', 'F3']),
        ('Lipos 2.7 nM',  ['E4', 'F4']),
        ('Lipos 1.4 nM',  ['E5', 'F5']),
        ('Lipos 0.68 nM',  ['E6', 'F6']),
        ('Lipos 0.34 nM',  ['E7', 'F7']),
        ('Lipos 0.17 nM',  ['E8', 'F8']),
        ('Lipos 0.085 nM',  ['E9', 'F9']),
        ('Lipos 0.042 nM',  ['E10', 'F10']),
        ('Lipos 0.021 nM',  ['E11', 'F11']),
        ('Lipos 0 nM',  ['E12', 'F12']),
        ('ANTS Lipos 12.1 nM, Bid/Bax', ['G1', 'G2']),
        ('ANTS Lipos 12.1 nM', ['G3', 'G4']),
        ])

lipo_conc_names = [
        'Lipos 0 nM',
        'Lipos 0.021 nM',
        'Lipos 0.042 nM',
        'Lipos 0.085 nM',
        'Lipos 0.17 nM',
        'Lipos 0.34 nM',
        'Lipos 0.68 nM',
        'Lipos 1.4 nM',
        'Lipos 2.7 nM',
        'Lipos 5.4 nM', 
        'Lipos 10.9 nM',
        'Lipos 21.7 nM',]
lipo_concs = [lc_name.split(' ')[1] for lc_name in lipo_conc_names]

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
preinc_file = os.path.abspath(os.path.join(data_path,
                                        '141016_Bax_depletion_preincubation.txt'))
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '141016_Bax_depletion_timecourse.txt'))
triton_file = os.path.abspath(os.path.join(data_path,
                                        '141016_Bax_depletion_triton.txt'))
bg_file = os.path.abspath(os.path.join(data_path,
                                     '141016_Bax_depletion_added_ANTS_EF.txt'))
# Assemble all the wells included in the layout
# http://stackoverflow.com/questions/406121/
# flattening-a-shallow-list-in-python
wells_to_read = list(itertools.chain(*layout.values()))

# Preincubation wells
preinc_wells = read_flexstation_kinetics(preinc_file)

# ANTS+liposomes to serve as background ctrls
bg_wells = read_flexstation_kinetics(bg_file)
bg_avgs = get_repeat_averages_by_well(bg_wells)

# Timecourse wells
timecourse_wells = read_flexstation_kinetics(timecourse_file)
"""The raw (unnormalized) timecourses."""
timecourse_wells = extract(wells_to_read, timecourse_wells)

# Initial timepoints
initial_wells_tc = get_first_points_by_well(timecourse_wells)

# Post-Triton values
final_wells = read_flexstation_kinetics(triton_file)
#final_wells = extract(wells_to_read, final_wells)
final_well_avgs = get_repeat_averages_by_well(final_wells)

# Averages of raw timecourses across replicates
(timecourse_averages, timecourse_stds) = averages(timecourse_wells, layout)
"""Averages of raw timecourses."""

# Normalized timecourses
norm_wells = get_normalized_well_timecourses(
        timecourse_wells, initial_wells_tc, final_well_avgs)
"""Timecourses normalized to min/max (initial/Triton) values."""

# Normalized timecourses, averaged
(norm_averages, norm_stds) = averages(norm_wells, layout)
"""Timecourses normalized and then averaged."""

# Get background average
#background = norm_averages['Bax 0 nM'][VALUE]

# Normalized and background subtracted
#bgsub_norm_wells = subtract_background(norm_wells, background)

# Normalized, background subtracted, averaged
#(bgsub_norm_averages, bgsub_norm_stds) = averages(bgsub_norm_wells, layout)

# First timepoint shifted to 0 (better for fitting)
#reset_bgsub_means = reset_first_timepoint_to_zero(bgsub_norm_averages)
"""Timecourses normalized, BG-subtracted, averaged, then with first point
shifted to t = 0."""
#reset_bgsub_sds = reset_first_timepoint_to_zero(bgsub_norm_stds)

# Pandas dataframe
#df = to_dataframe(bgsub_norm_averages, bgsub_norm_stds)
"""Pandas DataFrame version of reset_norm_means/sds"""

# Pore timecourses
#pores = get_average_pore_timecourses(bgsub_norm_averages)
"""Average pores derived by taking -log(1 - data)."""
def plot_preinc():
    # Get the baseline (untreated wells) and the treated wells
    ctrl_timecourses_all = extract(['G1', 'G2', 'G3', 'G4'], preinc_wells)
    # Plot the unnormalized timecourses along with baseline timecourses
    figure()
    plot_all(ctrl_timecourses_all)
    title('70 uL ANTS lipos, 100 nM Bid/Bax, preinc step')
    ylabel('RFU')
    xlabel('Time (sec)')

    # Normalized timecourses
    # Get baseline value from untreated ANTS wells
    baseline_wells = extract(['G3', 'G4'], preinc_wells)
    baseline_value = get_baseline_value(baseline_wells, num_timepoints=10)
    # Get just the treated ANTS wells
    ctrl_well_names = ['G1', 'G2']
    ctrl_timecourses = extract(ctrl_well_names, preinc_wells)
    ctrl_triton = extract(ctrl_well_names, final_well_avgs)
    # Normalize
    ctrl_norm_wells = get_normalized_well_timecourses(
            ctrl_timecourses, baseline_value, ctrl_triton)
    # Plot the normalized preinc timecourses
    plt.figure()
    plot_all(ctrl_norm_wells)
    plt.ylim([0, 1])
    plt.title('70 uL ANTS lipos, 100 nM Bid/Bax, preinc step')

    # Now take the average
    #ctrl_layout = collections.OrderedDict(
    #        [('Ctrls', ['G1', 'G2'])])
    #(ctrl_avgs, ctrl_stds) = averages(ctrl_norm_wells, ctrl_layout)
    #plt.errorbar(ctrl_avgs['Ctrls'][TIME], ctrl_avgs['Ctrls'][VALUE],
    #             yerr=ctrl_stds['Ctrls'][VALUE])
    #(ctrl_avgs2, ctrl_stds2) = averages(ctrl_norm_wells2, ctrl_layout)
    #ctrl_avgs2['Ctrls'][TIME] += 8000
    #ctrl_stds2['Ctrls'][TIME] += 8000
    #plt.errorbar(ctrl_avgs2['Ctrls'][TIME], ctrl_avgs2['Ctrls'][VALUE],
    #             yerr=ctrl_stds2['Ctrls'][VALUE])

    # Get the value of the ANTS wells during the main timecourse
    ctrl_timecourses2 = extract(ctrl_well_names, timecourse_wells)
    ctrl_norm_wells2 = get_normalized_well_timecourses(
                ctrl_timecourses2, baseline_value, ctrl_triton)

    # Plot the preinc and treatment timecourses on same plot
    plt.figure()
    for ctrl_well in ctrl_well_names:
        plt.plot(ctrl_norm_wells[ctrl_well][TIME],
                 ctrl_norm_wells[ctrl_well][VALUE], color='b')
        plt.plot(ctrl_norm_wells2[ctrl_well][TIME] + 8000,
                 ctrl_norm_wells2[ctrl_well][VALUE], color='g')
    plt.ylim([0, 1])
    plt.xlim([0, 2e4])
    plt.title('70 uL ANTS lipos, 100 nM Bid/Bax, step')

def plot_release_comparisons(plot_abs=True, plot_norm=True, bar_plot=True):
    preinc_fmaxes = []
    ref_fmaxes = []
    for lipo_conc_name in lipo_conc_names:
        # We'll take care of the calculations before plotting:
        ref_name = lipo_conc_name
        preinc_name = lipo_conc_name + ', Pre-inc'
        # Names of the wells
        preinc_well_names = layout[preinc_name]
        ref_well_names = layout[ref_name]
        # Get the averages for the ref wells from the background read
        bg_well_avgs = extract(ref_well_names, bg_avgs)
        # Take the mean of the two ref wells (rows E and F)
        bg_value = np.mean(bg_well_avgs.values())
        # Now that we have the background to subtract, normalize the wells
        # Normalized timecourses
        tc_wells_for_lipo_conc = extract(preinc_well_names + ref_well_names,
                                      timecourse_wells)
        final_wells_for_lipo_conc = extract(preinc_well_names + ref_well_names,
                                      final_well_avgs)
        norm_wells_for_lipo_conc = get_normalized_well_timecourses(
                tc_wells_for_lipo_conc, bg_value, final_wells_for_lipo_conc)

        # Now, to the plots.
        # Simplest version: plot the preinc and ref wells side by side, unnormed
        if plot_abs:
            plt.figure()
            for preinc_well in preinc_well_names:
                plt.plot(timecourse_wells[preinc_well][TIME],
                         timecourse_wells[preinc_well][VALUE], color='r')
            for ref_well in ref_well_names:
                plt.plot(timecourse_wells[ref_well][TIME],
                         timecourse_wells[ref_well][VALUE], color='b')
            # Plot the background value as a reference point on the unnormalized
            # plot.  This will help to determine if the background is wildly off
            # from the preincubated wells (since the background wells were taken
            # from the ref wells (rows E-F), not the preincubated wells, (rows
            # C-D))
            plt.hlines(bg_value, 0, 1e4, color='k', linestyle='--')
        # Plot the normalized version, where the baseline value from rows E and
        # F are used as the min value, and the per-well Triton reads are used
        # as the max value.
        if plot_norm:
            plt.figure()
            for preinc_well in preinc_well_names:
                plt.plot(norm_wells_for_lipo_conc[preinc_well][TIME],
                         norm_wells_for_lipo_conc[preinc_well][VALUE],
                         color='r')
            for ref_well in ref_well_names:
                plt.plot(norm_wells_for_lipo_conc[ref_well][TIME],
                         norm_wells_for_lipo_conc[ref_well][VALUE], color='b')
            plt.hlines(0, 0, 1e4, color='k', linestyle='--')
            plt.ylim([-0.1, 1.1])
            plt.xlim([-100, 1e4])
            plt.title(ref_name)

        # Ignore the 0.021 nM condition (because it's just like 0)
        for preinc_well in preinc_well_names:
            if preinc_well == 'C11' or preinc_well == 'D11':
                continue
            preinc_fmaxes.append(norm_wells_for_lipo_conc[preinc_well]
                                                         [VALUE][-1])
        for ref_well in ref_well_names:
            if ref_well == 'E11' or ref_well == 'F11':
                continue
            ref_fmaxes.append(norm_wells_for_lipo_conc[ref_well][VALUE][-1])

        # Make a pub-style plot of the 0.68 nM liposome condition
        if lipo_conc_name == 'Lipos 0.68 nM':
            set_fig_params_for_publication()
            plt.figure(figsize=(1.5, 1.5), dpi=300)
            for preinc_well in preinc_well_names:
                line_preinc = plt.plot(
                         norm_wells_for_lipo_conc[preinc_well][TIME],
                         norm_wells_for_lipo_conc[preinc_well][VALUE],
                         color=fig_purple)
            for ref_well in ref_well_names:
                line_ref = plt.plot(
                         norm_wells_for_lipo_conc[ref_well][TIME],
                         norm_wells_for_lipo_conc[ref_well][VALUE],
                         color=fig_orange)
            plt.subplots_adjust(left=0.25, bottom=0.23)
            ax = plt.gca()
            format_axis(ax)
            ax.set_xticks(np.linspace(0, 1e4, 6))
            ax.set_xticklabels([int(f) for f in np.linspace(0, 10, 6)])
            plt.ylim([0, 1])
            plt.xlim([-100, 1e4])
            plt.xlabel(r'Time (sec $\times 10^3$)', fontsize=fontsize)
            plt.ylabel(r'\% ANTS release', fontsize=fontsize)
            leg = plt.legend((line_ref[0], line_preinc[0]),
                       ('No preinc.', 'Preinc.'),
                       loc='right',
                       prop={'size': fontsize})
            leg.draw_frame(False)

    # Make the bar plot that is probably what we would show in a publication
    set_fig_params_for_publication()
    plt.figure(figsize=(4, 1), dpi=300)
    plt.subplots_adjust(bottom=0.34, right=0.95, top=0.87)
    loffset = 0.15
    width=0.15
    xcoords = np.zeros(len(ref_fmaxes))
    xcoords[0::2] = np.arange(len(ref_fmaxes) / 2)
    xcoords[1::2] = np.arange(len(ref_fmaxes) / 2) + width
    xcoords += loffset
    if bar_plot:
        rects_ref = plt.bar(xcoords, ref_fmaxes, width=width,
                            color=fig_orange)
        rects_preinc = plt.bar(xcoords + 2*width, preinc_fmaxes, width=width,
                               color=fig_purple)
        ax = gca()
        ax.set_xticks(xcoords[0::2] + 2*width)
        # Delete the 0.021 nM condition (because it's just like 0)
        lipo_concs_copy = copy(lipo_concs)
        del lipo_concs_copy[1]
        ax.set_xticklabels(lipo_concs_copy)
        ax.set_yticks(np.linspace(0, 1, 6)) # Ticks at 0, 0.2, 0.4, ... 1
        plt.xlim([0, len(ref_fmaxes) / 2])
        plt.ylim([np.min(preinc_fmaxes) * 1.5, 1])
        plt.xlabel('[Unlabeled liposomes] (nM)', fontsize=fontsize)
        plt.ylabel(r'\% ANTS release', fontsize=fontsize)
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_tick_params(which='both', labelsize=fontsize, pad=2,
                length=2, width=0.5)
        ax.yaxis.set_tick_params(which='both', direction='out',
                        labelsize=fontsize, pad=0, length=2, width=0.5)
        ax.xaxis.labelpad = 2
        ax.yaxis.labelpad = 2
        leg = plt.legend((rects_ref[0], rects_preinc[0]),
                   ('No preinc.', 'Preinc.'),
                   loc='upper right',
                   prop={'size': fontsize})
        leg.draw_frame(False)
        plt.show()

if __name__ == '__main__':
    # First, plot the ANTS controls to show that reaction should be at eq
    ion()
    plot_release_comparisons(plot_norm=False, plot_abs=False)
    for i in plt.get_fignums():
        plt.figure(i)
        plt.savefig('fig_141016_%d.pdf' % i)

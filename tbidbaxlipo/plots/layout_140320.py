from tbidbaxlipo.util.plate_assay import *
import itertools
import pickle
import collections
import sys
import os
import tbidbaxlipo.data
from matplotlib import pyplot as plt
import numpy as np
from tbidbaxlipo.util import fitting
from tbidbaxlipo.plots.titration_fits import TwoExp
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

if __name__ == '__main__':
    plt.ion()
    plot_data()
    sys.exit()
    #plot_lipo_background(lipo_bg_wells, lipo_bg_layout)

    # Try fitting the high conc trajectory
    y = bgsub_wells['A10'][VALUE]
    t = bgsub_wells['A10'][TIME]
    b = fitting.Parameter(np.mean(y[0:2]))
    fmax = fitting.Parameter(25.)
    k1 = fitting.Parameter(np.log(2)/2000.)
    k2 = fitting.Parameter(1e-3)
    k_bleach = fitting.Parameter(1.17e-5)
    def exp_func(t):
        return (b() + fmax()*(1 - np.exp(-k1()*(1 - np.exp(-k2()*t))*t))) * \
                np.exp(-k_bleach()*t)

    # One-parameter exp
    #def exp_func(t):
    #    return (b() + fmax()*(1 - np.exp(-k1()*t))) * \
    #            np.exp(-k_bleach()*t)

    fitting.fit(exp_func, [fmax, k1, k2], y, t)

    plt.figure()
    plt.plot(t, y, marker='o', linestyle='')
    plt.plot(t, exp_func(t), color='r', linewidth='2')

    print "Fmax: %f" % ((fmax() + b()) / b())

    sys.exit()
    #(fmax_arr, k1_arr, k2_arr, conc_list) = plot_two_exp_fits()
    #plot_fmax_curve(fmax_arr, conc_list)

    from tbidbaxlipo.plots import titration_fits
    """
    fit = TwoExp()
    pore_df = to_dataframe(pores, bgsub_norm_stds)
    fit.plot_fits_from_dataframe(pore_df)
    p = fit.fit_from_dataframe(pore_df)


    # Multiply Fmax values by molar concentration of liposomes
    concs = np.array(pore_df.columns.values, dtype='float')[1:-1]
    concentration_of_pores = p[1][1:-1] * 0.775
    plt.figure()
    plt.plot(concs, concentration_of_pores)

    # Fit a straight line to the concentrations
    from tbidbaxlipo.util import fitting
    m = fitting.Parameter(0.025)
    b = fitting.Parameter(0)
    def linear(x):
        return m()*x + b()
    fitting.fit(linear, [m, b], concentration_of_pores, concs)
    plt.plot(concs, linear(concs))
    """

    df = to_dataframe(norm_averages, norm_stds)
    concs = np.array(df.columns.values, dtype='float')
    fit = titration_fits.TwoExp()

    # Get background average
    background_time = norm_averages['Bax 0 nM'][TIME]
    background = norm_averages['Bax 0 nM'][VALUE]

    bg_rate = fit.fit_timecourse(background_time, background)

    plt.ion()
    plt.plot(background_time, background)
    plt.plot(background_time, fit.fit_func(background_time, bg_rate))
    t = np.linspace(0, 100000, 1000)
    plt.plot(t, fit.fit_func(t, bg_rate))
    import pdb; pdb.set_trace()

    fit = titration_fits.TwoExp()
    #fit.plot_fits_from_dataframe(subset_df)
    #p = fit.fit_from_dataframe(subset_df)
    fit.plot_fits_from_dataframe(df)
    p = fit.fit_from_dataframe(df)

    # With fitting of bg
    #print "bg_rate %f" % bg_rate
    fit = titration_fits.TwoExpWithBackground(bg_rate)
    #fit.plot_fits_from_dataframe(subset_df)
    #p = fit.fit_from_dataframe(subset_df)
    fit.plot_fits_from_dataframe(df)
    p = fit.fit_from_dataframe(df)

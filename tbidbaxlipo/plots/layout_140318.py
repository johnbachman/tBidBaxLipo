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

layout = collections.OrderedDict([
        ('Bax 185 nM, Lipos 1 mg/ml',  ['A1']),
        ('Bax 185 nM, Lipos 0.5 mg/ml',  ['A2']),
        ('Bax 185 nM, Lipos 0.25 mg/ml',  ['A3']),
        ('Bax 185 nM, Lipos 0.125 mg/ml',  ['A4']),
        ('Bax 185 nM, Lipos 0.063 mg/ml',  ['A5']),
        ('Bax 185 nM, Lipos 0.031 mg/ml',  ['A6']),
        ('Bax 185 nM, Lipos 0.016 mg/ml',  ['A7']),
        ('Bax 185 nM, Lipos 0.008 mg/ml',  ['A8']),
        ('Bax 185 nM, Lipos 0.004 mg/ml',  ['A9']),
        ('Bax 185 nM, Lipos 0.002 mg/ml',  ['A10']),
        ('Bax 185 nM, Lipos 0.001 mg/ml',  ['A11']),
        ('Bax 185 nM, Lipos 0 mg/ml',  ['A12']),
        ('Lipos 1 mg/ml',  ['B1']),
        ('Lipos 0.5 mg/ml',  ['B2']),
        ('Lipos 0.25 mg/ml',  ['B3']),
        ('Lipos 0.125 mg/ml',  ['B4']),
        ('Lipos 0.063 mg/ml',  ['B5']),
        ('Lipos 0.031 mg/ml',  ['B6']),
        ('Lipos 0.016 mg/ml',  ['B7']),
        ('Lipos 0.008 mg/ml',  ['B8']),
        ('Lipos 0.004 mg/ml',  ['B9']),
        ('Lipos 0.002 mg/ml',  ['B10']),
        ('Lipos 0.001 mg/ml',  ['B11']),
        ('Lipos 0 mg/ml',  ['B12']),
    ])
data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '140318_NBD_Bax_BimBH3_lipo_titration.txt'))

# Assemble all the wells included in the layout
# http://stackoverflow.com/questions/406121/
# flattening-a-shallow-list-in-python
wells_to_read = list(itertools.chain(*layout.values()))

# Timecourse wells
timecourse_wells = read_flexstation_kinetics(timecourse_file)
"""The raw (unnormalized) timecourses."""
timecourse_wells = extract(wells_to_read, timecourse_wells)

lipo_bg_conditions = [
        'Lipos 1 mg/ml',
        'Lipos 0.5 mg/ml',
        'Lipos 0.25 mg/ml',
        'Lipos 0.125 mg/ml',
        'Lipos 0.063 mg/ml',
        'Lipos 0.031 mg/ml',
        'Lipos 0.016 mg/ml',
        'Lipos 0.008 mg/ml',
        'Lipos 0.004 mg/ml',
        'Lipos 0.002 mg/ml',
        'Lipos 0.001 mg/ml',
        'Lipos 0 mg/ml']
lipo_bg_layout = extract(lipo_bg_conditions, layout)
lipo_bg_wells = extract([layout[cond][0] for cond in lipo_bg_conditions],
                        timecourse_wells)

bax_lipo_conditions = [
        'Bax 185 nM, Lipos 1 mg/ml',
        'Bax 185 nM, Lipos 0.5 mg/ml',
        'Bax 185 nM, Lipos 0.25 mg/ml',
        'Bax 185 nM, Lipos 0.125 mg/ml',
        'Bax 185 nM, Lipos 0.063 mg/ml',
        'Bax 185 nM, Lipos 0.031 mg/ml',
        'Bax 185 nM, Lipos 0.016 mg/ml',
        'Bax 185 nM, Lipos 0.008 mg/ml',
        'Bax 185 nM, Lipos 0.004 mg/ml',
        'Bax 185 nM, Lipos 0.002 mg/ml',
        'Bax 185 nM, Lipos 0.001 mg/ml',
        'Bax 185 nM, Lipos 0 mg/ml', ]
bax_lipo_layout = extract(bax_lipo_conditions, layout)
bax_lipo_wells = extract([layout[cond][0] for cond in bax_lipo_conditions],
                         timecourse_wells)

# Normalized and background subtracted
bgsub_wells = subtract_background_set(bax_lipo_wells, lipo_bg_wells)

#bgsub_norm_wells = subtract_background(norm_wells, background)

# Normalized, background subtracted, averaged
#(bgsub_norm_averages, bgsub_norm_stds) = averages(bgsub_norm_wells, layout)

# First timepoint shifted to 0 (better for fitting)
#reset_bgsub_means = reset_first_timepoint_to_zero(bgsub_norm_averages)
"""Timecourses normalized, BG-subtracted, averaged, then with first point
shifted to t = 0."""
#reset_bgsub_sds = reset_first_timepoint_to_zero(bgsub_norm_stds)


def plot_data():
    """Plots the data and various transformations of it."""
    ion()

    # Timecourse wells
    figure()
    plot_all(timecourse_wells)
    title("Raw timecourses")

    # Lipo bg wells
    figure()
    plot_all(lipo_bg_wells)
    title("Lipo background wells")

    # Background-subtracted wells
    figure()
    plot_all(bgsub_wells)
    title("Background-subtracted wells")

def plot_lipo_background(wells, layout):
    """Takes the lipo background well timecourses and plots the
    fluorescence values of the initial points as a function of lipo
    concentration."""
    init_vals = np.zeros(len(layout.keys()))
    conc_list = np.zeros(len(layout.keys()))
    for i, cond_name in enumerate(layout):
        conc = float(cond_name.split(' ')[1])
        well_name = layout[cond_name][0]
        init_val = wells[well_name][VALUE][0]
        init_vals[i] = init_val
        conc_list[i] = conc

    # Fit the liposome background to a line
    print "Fitting liposome background"
    m = fitting.Parameter(1.)
    b = fitting.Parameter(1.)
    def linear(s):
        return m()*s + b()
    result = fitting.fit(linear, [m, b], init_vals, conc_list)

    plt.figure()
    plt.plot(conc_list, init_vals, marker='o', linestyle='', color='b')
    plt.plot(conc_list, linear(conc_list), color='r')
    plt.xlabel('Liposome concentration (mg/ml)')
    plt.ylabel('RFU')
    plt.title('Liposome background fluorescence at t=0')
    print "m: %f" % m()
    print "b: %f" % b()
    return [m(), b()]

def plot_two_exp_fits():
    data = bgsub_wells
    fmax_arr = np.zeros((1, len(layout.keys())))
    k1_arr = np.zeros((1, len(layout.keys())))
    k2_arr = np.zeros((1, len(layout.keys())))
    conc_list = np.zeros(len(layout.keys()))

    for i, conc_str in enumerate(layout.keys()):
        plt.figure()
        wells = layout[conc_str]
        well_conc = float(conc_str.split(' ')[1])
        conc_list[i] = well_conc

        for j, well in enumerate(wells):
            well_data = data[well]
            time = np.array(well_data[TIME])
            y = np.array(well_data[VALUE])

            fit = TwoExp()
            (k1, fmax, k2) = fit.fit_timecourse(time, y)

            fmax_arr[j, i] = fmax
            k1_arr[j, i] = k1
            k2_arr[j, i] = k2

            plt.plot(time, y, 'b')
            plt.plot(time, fit.fit_func(time, (k1, fmax, k2)), 'r')
            plt.title(well_conc)
            plt.xticks([0, 2000, 4000, 6000, 8000])
            plt.ylabel('% Release')
            plt.xlabel('Time (sec)')

    plt.show()
    return (fmax_arr, k1_arr, k2_arr, conc_list)

def plot_k1_curve(k1_arr, conc_list):
    k1_means = np.mean(k1_arr, axis=0)
    k1_sds = np.std(k1_arr, axis=0)

    plt.figure()
    plt.errorbar(conc_list, k1_means, yerr= k1_sds / np.sqrt(3), color='r',
                 linestyle='', linewidth=2)
    plt.title('$k_1$')
    plt.xlabel('[Bax] (nM)')
    plt.ylabel('$k_1\ (\mathrm{sec}^{-1})$')
    """
    # Fit with exponential-linear curve
    vi = fitting.Parameter(0.05)
    vf = fitting.Parameter(0.025)
    tau = fitting.Parameter(0.02)
    # Define fitting function
    def biphasic(t):
        return (vf()*t) + ( (vi() - vf()) *
                            ((1 - np.exp(-tau()*t))/tau()) )
    fitting.fit(biphasic, [vi, vf, tau], k1_means, conc_list)
    plt.plot(conc_list, biphasic(conc_list), 'r', linewidth=2)

    # Fit with "double-binding curve"
    ksat = fitting.Parameter(0.0005)
    knonsat = fitting.Parameter(20)
    R0 = fitting.Parameter(20)
    def double_binding(t):
        return (R0() * t)/(ksat() + t) + (knonsat()*t)
    fitting.fit(double_binding, [R0, ksat, knonsat], k1_means, conc_list)
    plt.plot(conc_list, double_binding(conc_list), 'g', linewidth=2)

    plt.text(400, 0.00025,
            r'$\frac{R_0 Bax_0}{K_{sat} + Bax_0} + K_{nonsat} Bax_0$')
    plt.text(400, 0.00020, r'$R_0$ = %.3g nM' % R0())
    plt.text(400, 0.00015, r'$K_{sat}$ = %.3g nM' % ksat())
    plt.text(400, 0.00010, r'$K_{nonsat}$ = %.3g nM' % knonsat())

    plt.text(35, 0.00054,
            r'$v_f + (v_i - v_f) \left(\frac{1 - e^{-\tau t}}{\tau}\right)$')
    plt.text(35, 0.00049, r'$v_i$ = %.3g sec$^{-1}$ nM$^{-1}$' % vi())
    plt.text(35, 0.00044, r'$v_f$ = %.3g sec$^{-1}$ nM$^{-1}$' % vf())
    plt.text(35, 0.00039, r'$\frac{1}{\tau}$ = %.2f nM' % (1/tau()))
    plt.title('Biphasic fit of $k_1$')
    plt.show()
    """

def plot_k2_curve(k2_arr, conc_list):
    k2_means = np.mean(k2_arr, axis=0)
    k2_sds = np.std(k2_arr, axis=0)
    # Plot k2 on a linear scale
    plt.figure()
    plt.errorbar(conc_list, k2_means, yerr=k2_sds / np.sqrt(3))
    plt.title('$k_2$')

def plot_fmax_curve(fmax_arr, conc_list):
    fmax_means = np.mean(fmax_arr, axis=0)
    fmax_sds = np.std(fmax_arr, axis=0)

    # Plot of Fmax
    plt.figure()
    plt.ylabel(r'$F_{max}$ value (% Release)')
    plt.xlabel('[Bax] (nM)')
    plt.title(r'$F_{max}$')

    # Try fitting the fmax values to a hill function
    hill_vmax = fitting.Parameter(1)
    kd = fitting.Parameter(100)
    def hill(s):
        return ((hill_vmax() * s) / (kd() + s))
    fitting.fit(hill, [kd], fmax_means[1:], conc_list[1:])
    plt.plot(conc_list, hill(conc_list), 'r', linewidth=2, label='Hill fit')

    plt.text(35, 0.93, '$K_D$ = %.2f' % kd(), fontsize=16)
    #plt.text(35, 0.99, '$V_{max}$ = %.2f' % hill_vmax(), fontsize=16)

    # Try fitting the fmax values to an exponential function
    exp_fmax = fitting.Parameter(1.0)
    exp_k = fitting.Parameter(0.01)
    def exp_func(s):
        return exp_fmax() * (1 - np.exp(-exp_k()*s))
    fitting.fit(exp_func, [exp_fmax, exp_k], fmax_means[1:], conc_list[1:])
    plt.plot(conc_list, exp_func(conc_list), 'g', linewidth=2,
             label='Exp fit')

    # Plot the data
    plt.errorbar(conc_list[1:], fmax_means[1:],
                 yerr=fmax_sds[1:] / np.sqrt(3),
                 label='Data', linestyle='', linewidth=2)

    plt.legend(loc='lower right')
    plt.show()

if __name__ == '__main__':
    plt.ion()
    #plot_data()
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
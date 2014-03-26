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
from tbidbaxlipo.plots import titration_fits

def extract(keys, dict):
    extracted_dict = collections.OrderedDict([])
    for key in keys:
        extracted_dict[key] = dict[key]
    return extracted_dict


layout = collections.OrderedDict([
        ('Bax 0 nM',  ['A12', 'B12', 'C12']),
        ('Bax 0.8 nM', ['A11', 'B11', 'C11']),
        ('Bax 1.5 nM', ['A10', 'B10', 'C10']),
        ('Bax 3.1 nM', ['A09', 'B09', 'C09']),
        ('Bax 6.2 nM',   ['A08', 'B08', 'C08']),
        ('Bax 12.3 nM',  ['A07', 'B07', 'C07']),
        ('Bax 24.6 nM',  ['A06', 'B06', 'C06']),
        ('Bax 49.2 nM',  ['A05', 'B05', 'C05']),
        ('Bax 98.4 nM',  ['A04', 'B04', 'C04']),
        ('Bax 196.9 nM', ['A03', 'B03', 'C03']),
        ('Bax 393.8 nM', ['A02', 'B02', 'C02']),
        ('Bax 787.5 nM', ['A01', 'B01', 'C01']),
        ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '131210_Bax43C_PC_Timecourse.csv'))
#triton_file = os.path.abspath(os.path.join(data_path,
#                                        '130924_Bax_43C_Triton.csv'))

# Assemble all the wells included in the layout
# http://stackoverflow.com/questions/406121/
# flattening-a-shallow-list-in-python
wells_to_read = list(itertools.chain(*layout.values()))

# Timecourse wells
timecourse_wells = read_wallac(timecourse_file)
"""The raw (unnormalized) timecourses."""
timecourse_wells = extract(wells_to_read, timecourse_wells)

# Initial timepoints
initial_wells_tc = get_first_points_by_well(timecourse_wells)

# Post-Triton values
#final_wells = read_wallac(triton_file)
#final_wells = extract(wells_to_read, final_wells)
#final_well_avgs = get_repeat_averages_by_well(final_wells)

# Averages of raw timecourses across replicates
(timecourse_averages, timecourse_stds) = averages(timecourse_wells, layout)
"""Averages of raw timecourses."""

# Normalized timecourses
#norm_wells = get_normalized_well_timecourses(
#        timecourse_wells, initial_wells_tc, final_well_avgs)
"""Timecourses normalized to min/max (initial/Triton) values."""

# Normalized timecourses, averaged
#(norm_averages, norm_stds) = averages(norm_wells, layout)
"""Timecourses normalized and then averaged."""

# Get background average
background_time = timecourse_averages['Bax 0 nM'][TIME]
background = timecourse_averages['Bax 0 nM'][VALUE]

# Normalized and background subtracted
#bgsub_norm_wells = subtract_background(norm_wells, background)
bgsub_wells = subtract_background(timecourse_wells, background)

# Normalized, background subtracted, averaged
#(bgsub_norm_averages, bgsub_norm_stds) = averages(bgsub_norm_wells, layout)
(bgsub_averages, bgsub_stds) = averages(bgsub_wells, layout)

# First timepoint shifted to 0 (better for fitting)
#reset_bgsub_means = reset_first_timepoint_to_zero(bgsub_norm_averages)
"""Timecourses normalized, BG-subtracted, averaged, then with first point
shifted to t = 0."""
#reset_bgsub_sds = reset_first_timepoint_to_zero(bgsub_norm_stds)

# Pandas dataframe
#df = to_dataframe(bgsub_norm_averages, bgsub_norm_stds)
df = to_dataframe(bgsub_averages, bgsub_stds)
"""Pandas DataFrame version of reset_norm_means/sds"""

# Pore timecourses
#pores = get_average_pore_timecourses(bgsub_averages)
"""Average pores derived by taking -log(1 - data)."""


def plot_data():
    """Plots the data and various transformations of it."""
    ion()

    # Timecourse wells
    figure()
    plot_all(timecourse_wells)
    title("Raw timecourses")

    # Averages of raw timecourses across replicates
    figure()
    plot_all(timecourse_averages, errors=timecourse_stds)
    title("Raw timecourses, averaged")

    # Normalized timecourses
    #figure()
    #plot_all(norm_wells)
    #title("Normalized timecourses")

    # Normalized timecourses, background-subtracted
    figure()
    plot_all(bgsub_wells)
    #plot_all(bgsub_norm_wells)
    title("Normalized, BG-subtracted timecourses")

    # Normalized timecourses, averaged
    #figure()
    #plot_all(norm_averages, errors=norm_stds)
    #title("Normalized timecourses, averaged")

    # Normalized timecourses, background subtracted, averaged
    #figure()
    #plot_all(bgsub_norm_averages, errors=norm_stds)
    #title("Normalized timecourses, BG-subtracted, averaged")

    # First timepoint shifted to 0 (better for fitting)
    #figure()
    #plot_all(reset_bgsub_means)
    #title("Norm., BG-sub, avg., Reset to t = 0")

    # Pore timecourses
    #figure()
    #plot_all(pores)
    #title("Avg. pores per liposome")

def plot_two_exp_fits():
    data = bgsub_wells
    fmax_arr = np.zeros((3, len(layout.keys())))
    k1_arr = np.zeros((3, len(layout.keys())))
    k2_arr = np.zeros((3, len(layout.keys())))
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

            fit = titration_fits.TwoExp()
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
    #plot_two_exp_fits()

    from tbidbaxlipo.plots.titration_fits import OneExpFmax, TwoExp
    fit = OneExpFmax()
    #fit = TwoExp()
    fit.plot_fits_from_dataframe(df)
    p = fit.fit_from_dataframe(df)
    concs = fit.concs
    plt.figure()
    plt.plot(concs, p[0], marker='o', color='r')
    plt.title('k1')

    plt.figure()
    plt.plot(concs, p[1], marker='o', color='r')
    plt.title('Fmax')

    #plt.figure()
    #plt.plot(concs, p[2], marker='o', color='r')
    #plt.title('k2')
    sys.exit()

    # Pores
    #subset_pores = extract(active_concs, pores)
    #subset_pores_df = to_dataframe(subset_pores, subset_sds)

    #plot_all(subset_means)
    #plot_all(subset_pores)

    active_concs = ['Bax 7.9 nM',
                    'Bax 15.8 nM',
                    'Bax 31.5 nM',
                    'Bax 63 nM',
                    'Bax 126 nM']
    #subset_means = extract(active_concs, bgsub_norm_averages)
    #subset_means = extract(active_concs, norm_averages)
    #subset_sds = extract(active_concs, norm_stds)
    #subset_df = to_dataframe(subset_means, subset_sds)
    #concs = np.array(subset_df.columns.values, dtype='float')
    df = to_dataframe(norm_averages, norm_stds)
    concs = np.array(df.columns.values, dtype='float')

    fit = titration_fits.TwoExp()
    bg_rate = fit.fit_timecourse(background_time, background)
    print bg_rate
    plt.ion()
    plt.plot(background_time, background)
    plt.plot(background_time, fit.fit_func(background_time, bg_rate))
    t = np.linspace(0, 100000, 1000)
    plt.plot(t, fit.fit_func(t, bg_rate))

    fit = titration_fits.TwoExp()
    #fit.plot_fits_from_dataframe(subset_df)
    #p = fit.fit_from_dataframe(subset_df)
    fit.plot_fits_from_dataframe(df)
    p = fit.fit_from_dataframe(df)

    """
    plt.figure()
    plt.plot(concs, p[0], marker='o', color='r')
    plt.title('k1')

    plt.figure()
    plt.plot(concs, p[1], marker='o', color='r')
    plt.title('Fmax')

    plt.figure()
    plt.plot(concs, p[2], marker='o', color='r')
    plt.title('k2')
    """

    # With fitting of bg
    #print "bg_rate %f" % bg_rate
    fit = titration_fits.TwoExpWithBackground(bg_rate)
    #fit.plot_fits_from_dataframe(subset_df)
    #p = fit.fit_from_dataframe(subset_df)
    fit.plot_fits_from_dataframe(df)
    p = fit.fit_from_dataframe(df)

    """
    plt.figure()
    plt.plot(concs, p[0], marker='o', color='r')
    plt.title('k1')

    plt.figure()
    plt.plot(concs, p[1], marker='o', color='r')
    plt.title('Fmax')

    plt.figure()
    plt.plot(concs, p[2], marker='o', color='r')
    plt.title('k2')
    """

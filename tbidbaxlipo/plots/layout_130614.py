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
from scipy import stats
from tbidbaxlipo.util import fitting, set_fig_params_for_publication, \
            format_axis, fontsize
from tbidbaxlipo.plots import titration_fits

def extract(keys, dict):
    extracted_dict = collections.OrderedDict([])
    for key in keys:
        extracted_dict[key] = dict[key]
    return extracted_dict

layout = collections.OrderedDict([
        ('Bax 0 nM',  ['A12', 'B12', 'C12']),
        ('Bax 11 nM', ['A11', 'B11', 'C11']),
        ('Bax 16 nM', ['A10', 'B10', 'C10']),
        ('Bax 25 nM', ['A09', 'B09', 'C09']),
        ('Bax 37 nM',   ['A08', 'B08', 'C08']),
        ('Bax 55 nM',  ['A07', 'B07', 'C07']),
        ('Bax 83 nM',  ['A06', 'B06', 'C06']),
        ('Bax 124 nM',  ['A05', 'B05', 'C05']),
        ('Bax 187 nM',  ['A04', 'B04', 'C04']),
        ('Bax 280 nM', ['A03', 'B03', 'C03']),
        ('Bax 420 nM', ['A02', 'B02', 'C02']),
        ('Bax 630 nM', ['A01', 'B01', 'C01']),
        ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '130614_Bax_43C_Timecourse.csv'))
triton_file = os.path.abspath(os.path.join(data_path,
                                        '130614_Bax_43C_Triton.csv'))

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
final_wells = read_wallac(triton_file)
final_wells = extract(wells_to_read, final_wells)
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
background = norm_averages['Bax 0 nM'][VALUE]

# Normalized and background subtracted
bgsub_norm_wells = subtract_background(norm_wells, background)

# Normalized, background subtracted, averaged
(bgsub_norm_averages, bgsub_norm_stds) = averages(bgsub_norm_wells, layout)

# First timepoint shifted to 0 (better for fitting)
reset_bgsub_means = reset_first_timepoint_to_zero(bgsub_norm_averages)
"""Timecourses normalized, BG-subtracted, averaged, then with first point
shifted to t = 0."""
reset_bgsub_sds = reset_first_timepoint_to_zero(bgsub_norm_stds)

# Pandas dataframe
df = to_dataframe(bgsub_norm_averages, bgsub_norm_stds)
"""Pandas DataFrame version of reset_norm_means/sds"""

# Pore timecourses
pores = get_average_pore_timecourses(norm_averages)
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
    figure()
    plot_all(norm_wells)
    title("Normalized timecourses")

    # Normalized timecourses, background-subtracted
    figure()
    plot_all(bgsub_norm_wells)
    title("Normalized, BG-subtracted timecourses")

    # Normalized timecourses, averaged
    figure()
    plot_all(norm_averages, errors=norm_stds)
    title("Normalized timecourses, averaged")

    # Normalized timecourses, background subtracted, averaged
    figure()
    plot_all(bgsub_norm_averages, errors=norm_stds)
    title("Normalized timecourses, BG-subtracted, averaged")

    # First timepoint shifted to 0 (better for fitting)
    figure()
    plot_all(reset_bgsub_means)
    title("Norm., BG-sub, avg., Reset to t = 0")

    # Pore timecourses
    figure()
    plot_all(pores)
    title("Avg. pores per liposome")

def plot_two_exp_fits():
    data = bgsub_norm_wells
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

def calc_pore_size_by_poisson():
    set_fig_params_for_publication()

    (fmax_arr, k1_arr, k2_arr, conc_list) = \
         titration_fits.plot_two_exp_fits(bgsub_norm_wells, layout, plot=False)

    conc_lipos = 2.8

    fmax_means = np.mean(fmax_arr, axis=0)
    fmax_stds = np.std(fmax_arr, axis=0)
    ratios = conc_list / conc_lipos
    isf_pore_sizes = []
    k_max = 20
    for i, ratio in enumerate(ratios):
        if i == 0:
            continue
        # Probability of permeabilization is
        # 1 - poisson.cdf(pore_size, ratio)
        pore_sizes = np.arange(1, k_max)
        probs = [1 - stats.poisson.cdf(pore_size, ratio)
                 for pore_size in pore_sizes]
        isf_pore_size = stats.poisson.isf(fmax_means[i], ratio) + 1
        isf_pore_sizes.append(isf_pore_size)

        #plt.figure()
        #plt.plot(pore_sizes, probs, marker='o')
        #plt.xlabel('Pore size')
        #plt.ylim([-0.05, 1.05])
        #plt.title('Ratio of %f, isf pore size: %f' % (ratio, isf_pore_size))
        #sd_lines = [fmax_means[i] + fmax_stds[i]*num_sds
        #            for num_sds in np.arange(-2, 3)]
        #plt.hlines(sd_lines, 0, k_max, color='r')

    # Plot min pore size vs. Bax, with linear fit
    fig = plt.figure(figsize=(1.5, 1.5), dpi=300)
    ax = fig.gca()
    ax.plot(ratios[1:], isf_pore_sizes, marker='o', markersize=2,
            linestyle='') # Skip the 0 Bax pt
    ax.set_xlabel('[Bax]/[Lipo]')
    ax.set_ylabel('Predicted pore size')
    ax.set_xscale('log')
    ax.set_yscale('log')
    lbound = 0.1
    ubound = 200
    ax.set_xlim(lbound, 1000)
    #ax.set_ylim(1, 200)
    lin_fit = stats.linregress(ratios[1:], isf_pore_sizes)
    slope = lin_fit[0]
    intercept = lin_fit[1]
    ax.plot(ratios[1:], slope*ratios[1:] + intercept, color='r')
    interp_range = np.linspace(lbound, ratios[1])
    ax.plot(interp_range, slope*interp_range + intercept, color='r',
            linestyle='--', dashes=(2, 2))
    ax.plot(interp_range, [slope*ratios[1] + intercept] * len(interp_range),
            color='r', linestyle='--', dashes=(2, 2))
    #ax.set_title('Slope %f, intercept %f' % (slope, intercept),
    #             fontsize=fontsize)
    fig.subplots_adjust(left=0.22, bottom=0.19)
    format_axis(ax)



if __name__ == '__main__':
    plt.ion()
    #(fmax_arr, k1_arr, k2_arr, conc_list) = plot_two_exp_fits()
    #plot_fmax_curve(fmax_arr, conc_list)

    calc_pore_size_by_poisson()
    sys.exit()

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

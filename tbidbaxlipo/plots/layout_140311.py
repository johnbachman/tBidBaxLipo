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
from tbidbaxlipo.plots.titration_fits import TwoExp, Schwarz, Linear, OneExpFmax
import pymc
from scipy import stats

layout = collections.OrderedDict([
#        ('Buffer',    ['D7', 'D8', 'D9', 'D10', 'D11', 'D12']),
        ('Bax 0 nM',  ['A12', 'B12', 'C12']), #'D1', 'D2', 'D3', 'D4', 'D5', 'D6']),
        ('Bax 10.2 nM', ['A11', 'B11', 'C11']),
        ('Bax 15.4 nM', ['A10', 'B10', 'C10']),
        ('Bax 23 nM', ['A9', 'B9', 'C9']),
        ('Bax 35.5 nM',   ['A8', 'B8', 'C8']),
        ('Bax 51.8 nM',  ['A7', 'B7', 'C7']),
        ('Bax 77.7 nM',  ['A6', 'B6', 'C6']),
        ('Bax 117 nM',  ['A5', 'B5', 'C5']),
        ('Bax 175 nM',  ['A4', 'B4', 'C4']),
        ('Bax 262 nM', ['A3', 'B3', 'C3']),
        ('Bax 393 nM', ['A2', 'B2', 'C2']),
        ('Bax 590 nM', ['A1', 'B1', 'C1']),
        ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '140311_Bax_43C_Timecourse.txt'))
triton_file = os.path.abspath(os.path.join(data_path,
                                        '140311_Bax_43C_Triton2.txt'))

# Assemble all the wells included in the layout
# http://stackoverflow.com/questions/406121/
# flattening-a-shallow-list-in-python
wells_to_read = list(itertools.chain(*layout.values()))

# Timecourse wells
timecourse_wells = read_flexstation_kinetics(timecourse_file)
"""The raw (unnormalized) timecourses."""
timecourse_wells = extract(wells_to_read, timecourse_wells)

# Initial timepoints
initial_wells_tc = get_first_points_by_well(timecourse_wells)

# Post-Triton values
final_wells = read_flexstation_kinetics(triton_file)
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
#df = to_dataframe(bgsub_norm_averages, bgsub_norm_stds)
"""Pandas DataFrame version of reset_norm_means/sds"""

# Pore timecourses
#pores = get_average_pore_timecourses(bgsub_norm_averages)
pores = get_average_pore_timecourses(bgsub_norm_wells)
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

def plot_pore_fits(plot=True):
    data = bgsub_norm_wells
    k1_arr = np.zeros((3, len(layout.keys())))
    k2_arr = np.zeros((3, len(layout.keys())))
    tau_arr = np.zeros((3, len(layout.keys())))
    conc_list = np.zeros(len(layout.keys()))

    for i, conc_str in enumerate(layout.keys()):
        if plot:
            plt.figure()
        wells = layout[conc_str]
        well_conc = float(conc_str.split(' ')[1])
        conc_list[i] = well_conc

        for j, well in enumerate(wells):
            well_data = data[well]
            time = np.array(well_data[TIME])

            y = np.array(well_data[VALUE])
            y = -np.log(1 - y)

            fit = Schwarz()
            (k1, k2, tau) = fit.fit_timecourse(time, y)

            k1_arr[j, i] = k1
            k2_arr[j, i] = k2
            tau_arr[j, i] = tau

            if plot:
                plt.plot(time, y, 'b')
                plt.plot(time, fit.fit_func(time, (k1, k2, tau)), 'r')
                #plt.plot(time, fit.fit_func(time, (k1, fmax)), 'r')
                plt.title(well_conc)
                plt.xticks([0, 2000, 4000, 6000, 8000])
                plt.ylabel('% Release')
                plt.xlabel('Time (sec)')

    if plot:
        plt.show()
    return (k1_arr, k2_arr, tau_arr, conc_list)
    #return (fmax_arr, k1_arr, conc_list)

def plot_linear_fits(plot=True, max_time=50):
    data = bgsub_norm_wells
    k1_arr = np.zeros((3, len(layout.keys())))
    conc_list = np.zeros(len(layout.keys()))

    for i, conc_str in enumerate(layout.keys()):
        if plot:
            plt.figure()
        wells = layout[conc_str]
        well_conc = float(conc_str.split(' ')[1])
        conc_list[i] = well_conc

        for j, well in enumerate(wells):
            well_data = data[well]
            time = np.array(well_data[TIME][0:max_time])
            y = np.array(well_data[VALUE][0:max_time])
            fit = Linear()
            result = fit.fit_timecourse(time, y)
            k1 = result[0]

            k1_arr[j, i] = k1

            if plot:
                plt.plot(time, y, 'b')
                plt.plot(time, fit.fit_func(time, [k1]), 'r')
                plt.title(well_conc)
                plt.xticks([0, 2000, 4000, 6000, 8000])
                plt.ylabel('% Release')
                plt.xlabel('Time (sec)')

    if plot:
        plt.show()
    return (k1_arr, conc_list)

def plot_two_exp_fits(plot=True):
    data = pores
    #data = bgsub_norm_wells
    fmax_arr = np.zeros((3, len(layout.keys())))
    k1_arr = np.zeros((3, len(layout.keys())))
    k2_arr = np.zeros((3, len(layout.keys())))
    conc_list = np.zeros(len(layout.keys()))

    for i, conc_str in enumerate(layout.keys()):
        if plot:
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
            #fit = TwoExpFixedK2()
            #(k1, fmax) = fit.fit_timecourse(time, y)

            fmax_arr[j, i] = fmax
            k1_arr[j, i] = k1
            k2_arr[j, i] = k2

            if plot:
                plt.plot(time, y, 'b')
                plt.plot(time, fit.fit_func(time, (k1, fmax, k2)), 'r')
                #plt.plot(time, fit.fit_func(time, (k1, fmax)), 'r')
                plt.title(well_conc)
                plt.xticks([0, 2000, 4000, 6000, 8000])
                plt.ylabel('% Release')
                plt.xlabel('Time (sec)')

    if plot:
        plt.show()
    return (fmax_arr, k1_arr, k2_arr, conc_list)
    #return (fmax_arr, k1_arr, conc_list)

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

    # Try fitting the fmax values to a double-exponential
    exp_fmax1 = fitting.Parameter(1.)
    exp_fmax2 = fitting.Parameter(1.)
    exp_k1 = fitting.Parameter(0.001)
    exp_k2 = fitting.Parameter(0.001)
    def two_exp(s):
        return exp_fmax1()*(1 - np.exp(-exp_k1()*s)) + \
               exp_fmax2()*(1 - np.exp(-exp_k2()*s))
    fitting.fit(two_exp, [exp_fmax1, exp_fmax2, exp_k1, exp_k2], fmax_means[1:], conc_list[1:])
    plt.plot(conc_list, two_exp(conc_list), 'm', linewidth=2,
            label='Two-exp fit')

    # Try fitting the fmax values to a double-exponential
    hill_exp_fmax = fitting.Parameter(1.)
    hill_exp_kd = fitting.Parameter(100.)
    hill_exp_k = fitting.Parameter(0.01)
    def hill_exp(s):
        return hill_exp_fmax()*(1 - np.exp(
                                           ((-hill_exp_k() * s) / (hill_exp_kd() + s)) * s
                                          ))
    fitting.fit(hill_exp, [hill_exp_kd, hill_exp_k], fmax_means[1:], conc_list[1:])
    plt.plot(conc_list, hill_exp(conc_list), color='k', linewidth=2,
            label='Hill-exp fit')

    # Plot the data
    plt.errorbar(conc_list[1:], fmax_means[1:],
                 yerr=fmax_sds[1:] / np.sqrt(3),
                 label='Data', linestyle='', linewidth=2)

    plt.legend(loc='lower right')
    plt.show()

    # Pores per liposome at endpoint, hill fit
    plt.figure()
    import ipdb; ipdb.set_trace()
    pores_per_lipo = -np.log(1 - fmax_means[1:])
    plt.plot(conc_list[1:], pores_per_lipo, marker='o')
    plt.title('Pores per liposome at endpoint')

    hill_vmax = fitting.Parameter(1)
    kd = fitting.Parameter(100)
    def hill(s):
        return ((hill_vmax() * s) / (kd() + s))
    fitting.fit(hill, [kd, hill_vmax], pores_per_lipo, conc_list[1:])
    plt.plot(conc_list[1:], hill(conc_list[1:]), color='r')

def estimate_background_residuals():
    """Perform a fit of each background trajectory separately and then
    use these fits to estimate the variance of the residuals due to
    measurement error."""

    buffer_wells = extract(layout['Buffer'], timecourse_wells)
    plot_all(buffer_wells)

    (time_arr, buffer_reps) = get_replicates_for_condition(
                                    timecourse_wells, layout, 'Buffer')
    print "bg_init   k_decay   k_linear   bg_bg"
    num_timepoints = len(time_arr)
    num_replicates = buffer_reps.shape[1]
    residuals = np.zeros((num_timepoints, num_replicates))

    param_matrix = np.zeros((num_replicates, 4))

    for i in range(num_replicates):
        plt.figure()
        plt.plot(time_arr, buffer_reps[:,i])

        bg_init = fitting.Parameter(1.)
        k_decay = fitting.Parameter(1e-2)
        bg_bg = fitting.Parameter(5.6)
        k_linear = fitting.Parameter(1e-4)

        def exp_decay(t):
            return bg_init()*np.exp(-k_decay() * t) - k_linear()*t + bg_bg()

        (residuals[:,i], result) = fitting.fit(exp_decay, [bg_init, k_decay, bg_bg, k_linear],
                                               buffer_reps[:,i], time_arr)
        plt.plot(time_arr, exp_decay(time_arr), linewidth=2, color='r')
        print bg_init(), k_decay(), k_linear(), bg_bg()
        param_matrix[i,:] = [bg_init(), k_decay(), k_linear(), bg_bg()]

    plt.figure()
    plt.hist(residuals.reshape(num_replicates * num_timepoints))
    print "Mean"
    print np.mean(residuals)
    print "Variance"
    print np.var(residuals)
    print "SD"
    print np.std(residuals)
    print param_matrix
    print "Mean"
    print np.mean(param_matrix, axis=0)
    print "SD"
    print np.std(param_matrix, axis=0)

    return np.var(residuals)

def estimate_bg_release():
    lipo_wells = extract(layout['Bax 0 nM'], timecourse_wells)
    plot_all(lipo_wells)

    (time_arr, lipo_reps) = get_replicates_for_condition(
                                    timecourse_wells, layout, 'Bax 0 nM')
    num_timepoints = len(time_arr)
    num_replicates = lipo_reps.shape[1]
    residuals = np.zeros((num_timepoints, num_replicates))
    param_matrix = np.zeros((num_replicates, 2))

    avg_bg = timecourse_averages['Buffer'][VALUE]
    for i in range(num_replicates):
        plt.figure()
        plt.plot(time_arr, lipo_reps[:,i])
        bg_sub = lipo_reps[:,i] - avg_bg
        plt.plot(time_arr, bg_sub)

        lipo_init = fitting.Parameter(20.)
        k_app = fitting.Parameter(1e-2)

        def linear(t):
            return lipo_init() + k_app()*t

        (residuals[:,i], result) = fitting.fit(linear, [lipo_init, k_app],
                                               bg_sub, time_arr)
        plt.plot(time_arr, linear(time_arr), linewidth=2, color='r')
        print lipo_init(), k_app()
        param_matrix[i,:] = [lipo_init(), k_app()]

    plt.figure()
    plt.hist(residuals.reshape(num_replicates * num_timepoints))
    print "Mean"
    print np.mean(residuals)
    print "Variance"
    print np.var(residuals)
    print "SD"
    print np.std(residuals)
    print param_matrix
    print "Mean"
    print np.mean(param_matrix, axis=0)
    print "SD"
    print np.std(param_matrix, axis=0)

    return

def approach_2(avgs, stds, layout):
    """In this approach we perform MCMC fitting of the replicate averages
    for each concentration and plot the uncertainty associated with each
    on the final titration plot."""
    # Get the average trajectories, with standard errors
    #(avgs, stderrs) = averages(wells, layout, stderr=True)

    iter=50000
    burn=iter/4.
    thin=50.
    num_samples = (iter - burn) / thin

    num_concs = len(layout.keys())
    k1_arr = np.zeros((num_concs, 2))
    k2_arr = np.zeros((num_concs, 2))
    fmax_arr = np.zeros((num_concs, 2))
    conc_list = np.zeros(num_concs)
    k1_traces = np.zeros((num_concs, num_samples))

    def exp_func(k1_, k2_, fmax_, t):
        """Two-exponential fitting function."""
        return (fmax_ * (1 - np.exp(-k1_ *
                                      (1 - np.exp(-k2_ * t)) * t)))

    plt.figure()
    num_rows = 4
    num_cols = 3
    for i, conc_str in enumerate(layout.keys()):
        print conc_str
        conc = float(conc_str.split(' ')[1])
        conc_list[i] = conc
        time = avgs[conc_str][TIME]
        curve = avgs[conc_str][VALUE]
        plt.subplot(num_rows, num_cols, i+1)
        plt.errorbar(time, curve, yerr=stds[conc_str][VALUE], color='r')

        k1 = pymc.Uniform('k1', 0, 1 / 24.)
        k2 = pymc.Uniform('k2', 0, 1 / 24.)
        fmax = pymc.Uniform('fmax', 0, 1)

        @pymc.stochastic(plot=True, observed=True, verbose=0)
        def exp_fit(k1_=k1, k2_=k2, fmax_=fmax, value=curve):
            f = exp_func(k1_, k2_, fmax_, time)
            #cv = 0.01 # Magic number, should estimate from data
            #variance = (np.min(curve) * cv)**2
            variance = stds[conc_str][VALUE] **2
            # a hack b/c otherwise the variance in the normalized data is 0
            variance[0] = 0.014
            log_lkl = -np.sum(((curve - f)**2) / (2 * variance))
            return log_lkl

        mcmc = pymc.MCMC([k1, k2, fmax, exp_fit])
        mcmc.sample(iter=iter, burn=burn, thin=thin)
        print

        # TODO: should select randomly, not just last n samples
        num_to_plot = 500
        k1_trace = mcmc.trace('k1')[:]
        k2_trace = mcmc.trace('k2')[:]
        fmax_trace = mcmc.trace('fmax')[:]
        for j in range(num_to_plot):
            pass
            plt.plot(time, exp_func(k1_trace[-j], k2_trace[-j],
                                    fmax_trace[-j], time),
                     alpha=0.05, color='b')

        #pymc.Matplot.plot(mcmc)

        plt.title(conc_str)
        plt.xlabel('Time (sec)')
        plt.ylabel('Pct release')
        #import ipdb; ipdb.set_trace()
        k1_arr[i] = [np.mean(k1_trace), np.std(k1_trace)]
        k2_arr[i] = [np.mean(k2_trace), np.std(k2_trace)]
        fmax_arr[i] = [np.mean(fmax_trace), np.std(fmax_trace)]
        #min_rfu_arr[i] = [np.mean(min_trace), np.std(min_trace)]
        #k_arr[i] = [np.mean(k_trace), np.std(k_trace)]
        k1_traces[i, :] = k1_trace

    plt.tight_layout(pad=0, w_pad=0, h_pad=0)

    # Plot values of fitted parameters
    plt.figure()
    plt.errorbar(conc_list, k1_arr[:,0], yerr=k1_arr[:,1])
    plt.figure()
    plt.errorbar(conc_list, k2_arr[:,0], yerr=k2_arr[:,1])
    plt.figure()
    plt.errorbar(conc_list, fmax_arr[:,0], yerr=fmax_arr[:,1])
    plt.figure()
    plt.boxplot(k1_traces[:-1].T)

def plot_k_curves(k_arr, conc_list):
    k_means = np.mean(k_arr, axis=0)
    k_sds = np.std(k_arr, axis=0)

    plt.figure()
    plt.errorbar(conc_list, k_means, yerr= k_sds / np.sqrt(3), color='r',
                 linestyle='', linewidth=2)
    plt.title('$k$')
    plt.xlabel('[Bax] (nM)')
    plt.ylabel('$k\ (\mathrm{sec}^{-1})$')

def plot_initial_slope_curves(k_arr, conc_list):
    k_means = np.mean(k_arr, axis=0)
    k_sds = np.std(k_arr, axis=0)

    plt.figure()
    plt.errorbar(conc_list, k_means, yerr= k_sds / np.sqrt(3), color='r',
                 linestyle='', linewidth=2)
    plt.title('$k$')
    plt.xlabel('[Bax] (nM)')
    plt.ylabel('$k\ (\mathrm{sec}^{-1})$')

    fmax = fitting.Parameter(1e-3)
    k = fitting.Parameter(np.log(2)/300.)
    def exp_func(s):
        return fmax()*(1 - np.exp(-k()*s))
    fitting.fit(exp_func, [fmax, k], k_means, conc_list)
    plt.plot(conc_list, exp_func(conc_list))

    k_linear = fitting.Parameter(4e-7)
    def linear(s):
        return k_linear()*s
    fitting.fit(linear, [k_linear], k_means[0:7], conc_list[0:7])
    plt.plot(conc_list, linear(conc_list), color='b')
    return k_linear()

if __name__ == '__main__':
    plt.ion()
    #approach_2(bgsub_norm_averages, bgsub_norm_stds, layout)
    plot_data()
    sys.exit()
    #(k1_arr, k2_arr, tau_arr, conc_list) = plot_pore_fits(plot=True)
    #plot_k_curves(k1_arr, conc_list)
    #plot_k_curves(k2_arr, conc_list)
    #plot_k_curves(tau_arr, conc_list)

    #(k1_arr, conc_list) = plot_linear_fits(plot=False)
    #bax_rate_slope = plot_initial_slope_curves(k1_arr, conc_list)

    (fmax_arr, k1_arr, k2_arr, conc_list) = plot_two_exp_fits(plot=False)
    plot_fmax_curve(fmax_arr, conc_list)
    plot_k1_curve(k1_arr, conc_list)
    plot_k2_curve(k2_arr, conc_list)

    fmax_means = np.mean(fmax_arr, axis=0)
    k1_means = np.mean(k1_arr, axis=0)
    k2_means = np.mean(k2_arr, axis=0)

    def two_exp_deriv(t, k1, fmax, k2):
        return np.exp(-(k1 - np.exp(-k2*t)*k1 + k2)*t) * \
                fmax * k1 * \
                (-1 + np.exp(k2*t) + k2*t)

    def one_exp_deriv(t, k1, fmax):
        return np.exp(-k1 *t)*fmax*k1

    # One Exp fits to 10.2 nM
    te = OneExpFmax()
    t = bgsub_norm_wells['A11'][TIME]
    y = bgsub_norm_wells['A11'][VALUE]
    params = te.fit_timecourse(t, y)

    plt.figure()
    plt.plot(t, y)
    plt.plot(t, te.fit_func(t, params))

    #plt.figure()
    #plt.plot(t, one_exp_deriv(t, *params))
    plt.figure()
    plt.plot(t, one_exp_deriv(t, *params) / bax_rate_slope)

    t_start = 0
    t_end = 100
    amt_bax_lost = one_exp_deriv(t_start, *params) / bax_rate_slope - \
                   one_exp_deriv(t_end, *params) / bax_rate_slope
    amt_permeabilized = te.fit_func(t_end, params) - \
                        te.fit_func(t_start, params)
    conc_lipos = 5.16
    k_lambda = amt_bax_lost / conc_lipos
    pore_size = stats.poisson.isf(amt_permeabilized, k_lambda)
    print pore_size
    sys.exit()

    plt.figure()
    plt.plot(t, one_exp_deriv(t, 1e-4, 4e-2))
    #plt.plot(t, two_exp_deriv(t, 1e-4, 4e-2, 1), color='r')

    fmax_means = np.mean(fmax_arr, axis=0)
    fmax_stds = np.std(fmax_arr, axis=0)
    ratios = conc_list / 5.16

    derivs = np.zeros(len(fmax_means)-1)
    for i in np.arange(1, len(fmax_means)):
        if i == 1:
            prev_mean = 0
        else:
            prev_mean = fmax_means[i-1]

        fmax_diff = fmax_means[i] - prev_mean
        conc_diff = conc_list[i] - conc_list[i-1]
        derivs[i-1] = fmax_diff / float(conc_diff)
    plt.figure()
    plt.plot(conc_list[1:], derivs, marker='o')

    exp_lambda = fitting.Parameter(1/20.)
    def exp_dist(x):
        return exp_lambda() * np.exp(-exp_lambda() * x)
    (residuals, result) = fitting.fit(exp_dist, [exp_lambda], derivs, conc_list[1:])
    plt.plot(conc_list[1:], exp_dist(conc_list[1:]))
    sys.exit()

    """
    k_max = 20
    for i, ratio in enumerate(ratios):
        if i == 0:
            continue
        # Probability of permeabilization is
        # 1 - poisson.cdf(pore_size, ratio)
        pore_sizes = np.arange(1, k_max)
        probs = [1 - stats.poisson.cdf(pore_size, ratio) for pore_size in pore_sizes]
        plt.figure()
        plt.plot(pore_sizes, probs)
        plt.xlabel('Pore size')
        plt.ylim([-0.05, 1.05])
        plt.title('Ratio of %f' % ratio)
        sd_lines = [fmax_means[i] + fmax_stds[i]*num_sds for num_sds in np.arange(-2, 3)]
        plt.hlines(sd_lines, 0, k_max, color='r')
    """

    # Predicted Fmax curve for fixed pore size
    pore_size = 1
    probs = [1 - stats.poisson.cdf(pore_size, ratio/10.) for ratio in ratios]
    plt.figure()
    plt.errorbar(conc_list[1:], fmax_means[1:], yerr=fmax_stds[1:])
    plt.plot(conc_list, probs)
    sys.exit()

    cov_matrix = result[1] * res_var

    bg_init_sd = np.sqrt(cov_matrix[0,0])
    k_decay_sd = np.sqrt(cov_matrix[1, 1])
    bg_bg_sd = np.sqrt(cov_matrix[2, 2])

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

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
import string
import pymc

import matplotlib
font = {'size': 8}
matplotlib.rc('font', ** font)

def extract(keys, dict):
    extracted_dict = collections.OrderedDict([])
    for key in keys:
        extracted_dict[key] = dict[key]
    return extracted_dict

layout = collections.OrderedDict([
        #('Liposomes 1.33 mg/ml', ['%s01' % string.uppercase[i] for i in range(8)]),
        #('Liposomes 0.89 mg/ml', ['%s02' % string.uppercase[i] for i in range(8)]),
        ('Liposomes 0.59 mg/ml', ['%s03' % string.uppercase[i] for i in range(8)]),
        ('Liposomes 0.40 mg/ml', ['%s04' % string.uppercase[i] for i in range(8)]),
        ('Liposomes 0.26 mg/ml', ['%s05' % string.uppercase[i] for i in range(8)]),
        ('Liposomes 0.18 mg/ml', ['%s06' % string.uppercase[i] for i in range(8)]),
        ('Liposomes 0.12 mg/ml', ['%s07' % string.uppercase[i] for i in range(8)]),
        ('Liposomes 0.078 mg/ml', ['%s08' % string.uppercase[i] for i in range(8)]),
        ('Liposomes 0.052 mg/ml', ['%s09' % string.uppercase[i] for i in range(8)]),
        ('Liposomes 0.035 mg/ml', ['%s10' % string.uppercase[i] for i in range(8)]),
        ('Liposomes 0.023 mg/ml', ['%s11' % string.uppercase[i] for i in range(8)]),
        ('Liposomes 0 mg/ml', ['%s12' % string.uppercase[i] for i in range(8)]),
        ])

data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                                        '140305_ANTS_bg_envision.csv'))

#timecourse_file = os.path.abspath(os.path.join(data_path,
#                                        '140220_ANTS_plt3881_buffer.csv'))
#timecourse_file = os.path.abspath(os.path.join(data_path,
#                                        '140220_ANTS_bg_release.csv'))

#triton_file = os.path.abspath(os.path.join(data_path,
#                                        '130614_Bax_43C_Triton.csv'))

# Assemble all the wells included in the layout
# http://stackoverflow.com/questions/406121/
# flattening-a-shallow-list-in-python
wells_to_read = list(itertools.chain(*layout.values()))

# Timecourse wells
timecourse_wells = read_envision(timecourse_file)
"""The raw (unnormalized) timecourses."""
timecourse_wells = extract(wells_to_read, timecourse_wells)

# Initial timepoints
#initial_wells_tc = get_first_points_by_well(timecourse_wells)

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
#background = timecourse_averages['Buffer'][VALUE]

# Normalized and background subtracted
#bgsub_wells = subtract_background(timecourse_averages, background)

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
#pores = get_average_pore_timecourses(norm_averages)
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

    return

    # Background-subtracted
    figure()
    plot_all(bgsub_wells)
    title("BG-subtracted timecourses")

    # Normalized timecourses
    figure()
    plot_all(norm_wells)
    title("Normalized timecourses")


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

class FitTest(object):

    def __init__(self, k_mean, k_sd, initial_rfu_mean, initial_rfu_cv,
            pct_increase_mean, pct_increase_sd, read_cv, num_timepoints):
        self.k_mean = k_mean
        self.k_sd = k_sd
        self.initial_rfu_mean = initial_rfu_mean
        self.initial_rfu_cv = initial_rfu_cv
        self.pct_increase_mean = pct_increase_mean
        self.pct_increase_sd = pct_increase_sd
        self.read_cv = read_cv
        self.num_timepoints = num_timepoints
        self.t = np.linspace(0, 12000, num_timepoints)

    def get_synthetic_replicates(self, num_replicates, plot=False):
        reps = np.zeros((num_replicates, self.num_timepoints))
        # Time vect
        if plot:
            plt.figure()

        for i in range(num_replicates):
            # The CV for the initial value (e.g., due to pipetting error seems
            # to be 1-2%, we'll say 1.3% (the mean of the CVs across concentration)
            # We'll use a mean of 400,000 RFUs so we're in a comparable range
            #initial_rfu_mean = 4e5
            #initial_rfu_cv = 0.013
            initial_rfu_sd = self.initial_rfu_mean * self.initial_rfu_cv
            initial_rfu = initial_rfu_sd * np.random.randn(1) + self.initial_rfu_mean
            # We'll assume a similar CV for the within replicate read error
            #read_cv = 0.01
            read_sd = (self.initial_rfu_mean * self.read_cv)
            # Mean and variance for the kinetic constant
            #k_mean = 1.5e-4
            #k_sd = 2.5e-5
            k = self.k_sd * np.random.randn(1) + self.k_mean
            # For the percent increase
            #pct_increase_mean = 0.1
            #pct_increase_sd = 0.03
            pct_increase = self.pct_increase_sd * np.random.randn(1) \
                           + self.pct_increase_mean
            # Generate trajectory
            trajectory = (pct_increase * initial_rfu)*(1 - np.exp(-k * self.t)) \
                         + initial_rfu
            noise_vector = read_sd * np.random.randn(self.num_timepoints)
            trajectory_with_noise = trajectory + noise_vector
            reps[i,:] = trajectory_with_noise

            if plot:
                plt.plot(self.t, trajectory_with_noise, color='b', alpha=0.2)
        avg_trajectory = np.mean(reps, axis=0)

        # Plot the average
        if plot:
            plt.errorbar(self.t, avg_trajectory, \
                         yerr=np.std(reps, axis=0)/np.sqrt(num_replicates),
                         color='r', linewidth=2)
        return reps

    def fit_average_by_leastsq(self, reps):
        """Fit the average with an exponential model

        Returns a dict with the parameter values for the fit.
        """

        avg_trajectory = np.mean(reps, axis=0)

        k_p = fitting.Parameter(np.log(2)/6000.)
        initial_rfu_p = fitting.Parameter(avg_trajectory[0])
        pct_increase_p = fitting.Parameter((avg_trajectory[-1]-avg_trajectory[0])/
                                            avg_trajectory[0])

        def exp_func(t):
            return (pct_increase_p() * initial_rfu_p()) * (1 - np.exp(-k_p() * t)) + \
                    initial_rfu_p()
        fitting.fit(exp_func, [k_p, initial_rfu_p, pct_increase_p], avg_trajectory, self.t)
        return {'k': k_p(),
                'initial_rfu': initial_rfu_p(),
                'pct_increase': pct_increase_p()}

    def leastsq_error(self, num_trials, num_replicates, bins=20):
        p_errs = np.zeros((num_trials, 3))
        for i in range(num_trials):
            reps = self.get_synthetic_replicates(num_replicates, plot=False)
            p_fits = self.fit_average(reps)
            p_errs[i,:] = [p_fits['k'], p_fits['initial_rfu'],
                           p_fits['pct_increase']]
        plt.figure()
        plt.subplot(1, 3, 1)
        plt.hist(p_errs[:,0], bins=bins)
        plt.title('k')
        plt.subplot(1, 3, 2)
        plt.hist(p_errs[:,1], bins=bins)
        plt.title('initial_rfu')
        plt.subplot(1, 3, 3)
        plt.hist(p_errs[:,2], bins=bins)
        plt.title('pct_increase')

    def fit_average_by_mcmc(self, reps, bins=20):
        iter=50000
        burn=iter/4.
        thin=5.
        num_samples = (iter - burn) / thin

        avg_trajectory = np.mean(reps, axis=0)
        stderr_trajectory = np.std(reps, axis=0) / np.sqrt(reps.shape[0])
        #num_concs = len(layout.keys())
        #min_rfu_arr = np.zeros((num_concs, 2))
        #max_rfu_arr = np.zeros((num_concs, 2))
        #k_arr = np.zeros((num_concs, 2))
        #conc_list = np.zeros(num_concs)
        #k_traces = np.zeros((num_concs, num_samples))

        def exp_func(pct, min, k, time):
            return pct * min * (1 - np.exp(-k * time)) + min

        # Priors
        pct_increase = pymc.Uniform('pct_increase', 0, 1)
        initial_rfu = pymc.Normal('initial_rfu', mu=avg_trajectory[0],
                              tau=1/(stderr_trajectory[0]**2))
        # Range between the min and max timescale
        # TODO Should get rid of magic numbers, calculate from data instead
        k = pymc.Uniform('k', 1 / 12000., 1 / 120.)

        @pymc.stochastic(plot=True, observed=True, verbose=0)
        def exp_fit(pct=pct_increase, min=initial_rfu, k=k, value=avg_trajectory):
            f = exp_func(pct, min, k, self.t)
            #cv = 0.01 # Magic number, should estimate from data
            #variance = (np.min(curve) * cv)**2
            variance = stderr_trajectory ** 2
            log_lkl = -np.sum(((avg_trajectory - f)**2) / (2 * variance))
            return log_lkl

        mcmc = pymc.MCMC([pct_increase, initial_rfu, k, exp_fit])
        mcmc.sample(iter=iter, burn=burn, thin=thin)
        print

        # TODO: should select randomly, not just last n samples
        plt.figure()
        plt.errorbar(self.t, avg_trajectory, yerr=stderr_trajectory, color='r')
        num_to_plot = 500
        pct_trace = mcmc.trace('pct_increase')[:]
        initial_rfu_trace = mcmc.trace('initial_rfu')[:]
        k_trace = mcmc.trace('k')[:]
        for j in range(num_to_plot):
            plt.plot(self.t,
                     exp_func(pct_trace[-j], initial_rfu_trace[-j], k_trace[-j],
                              self.t),
                     alpha=0.05, color='b')

        #pymc.Matplot.plot(mcmc)
        plt.figure()
        plt.subplot(1, 3, 1)
        plt.hist(k_trace, bins=bins)
        plt.title('k')
        plt.subplot(1, 3, 2)
        plt.hist(initial_rfu_trace, bins=bins)
        plt.title('initial_rfu')
        plt.subplot(1, 3, 3)
        plt.hist(pct_trace, bins=bins)
        plt.title('pct_increase')

    # I. Working with averages of replicates only

    # 1. Leastsq-fit of average.  Fitting the average of even as little as five
    # replicates does reasonably well in that it generally gets the true mean
    # of each of the parameters. However, there is error in these estimates and
    # since for any single experiment it produces only a single value we have
    # no measure of how inaccurate our estimate may be. This would be useful,
    # for example, in comparing the kinetic constants between different
    # concentrations. It may be that there is a theory, as there is for linear
    # regression, connecting the error in the data points themselves to the
    # errors in the parameters, but I don't know what that is or if it would
    # generalize to more complicated models.

    # 2. Running MCMC on average.
    # Running MCMC on an average of even 5 or 10 reps produces an estimate of the
    # uncertainty in the parameters. Naturally, the uncertainty goes down as the
    # number of replicates goes up (and hence the standard error goes down). This
    # allows for good estimates of the means underlying the various distributions.

    # So, is this adequate? Is it important to know the variance in the
    # replicate distribution or is it adequate to simply know the uncertainty
    # associated with the average?  I suppose the issue is that if you don't
    # know the spread associated with single replicates, you don't know the
    # range of possibilities for the background associated with a particular
    # well--you will simply subtract use the mean value, with any associated
    # uncertainty.  It doesn't make sense that the variation in background
    # release that may be affecting your experimental condition should be
    # affected by the number of replicates you do to estimate the mean. In
    # fact, for this reason it would seem that what you really want to do is
    # estimate the mean and variance associated with the individual wells.  For
    # this reason, it would seem desirable to work with fits of the original
    # replicates.

    def fit_replicates_by_leastsq(self, reps):
        # Iterate over the replicates
        # For each one, run a fit (inside inner loop)
        # Record values of parameters
        # At end, plot mean, variance, etc. 
        pass

    def fit_replicates_by_mcmc(self, reps):
        # Iterate over the replicates
        # For each one, run a fit
        # Record complete traces
        # Compare two approaches for calculating mean and variance
        #   - Take means of each fit, and calculate mean and var across
        #   - Calculate mean and variance by sampling from traces many times
        pass

param_dict = {'k_mean': 1.5e-4, 'k_sd': 2.5e-5,
              'initial_rfu_mean': 4e5, 'initial_rfu_cv': 0.013,
              'pct_increase_mean': 0.1, 'pct_increase_sd': 0.03,
              'read_cv': 0.01, 'num_timepoints': 101}

def approach_1(wells, layout):
    """In this approach, we take the mean of each condition and then perform
    a single fit to the mean. We then plot the value of each parameter
    as a function of the concentration to look for patterns.

    Note that the data is expected to be in raw format, not pre-averaged.
    """
    (avgs, stderrs) = averages(wells, layout, stderr=True)

    num_concs = len(layout.keys())
    min_rfu_arr = np.zeros((num_concs, 2))
    max_rfu_arr = np.zeros((num_concs, 2))
    k_arr = np.zeros((num_concs, 2))
    conc_list = np.zeros(num_concs)

    plt.figure()
    num_rows = 4
    num_cols = 3
    for i, conc_str in enumerate(layout.keys()):
        print conc_str
        conc = float(conc_str.split(' ')[1])
        conc_list[i] = conc
        time = avgs[conc_str][TIME][0:100]
        curve = avgs[conc_str][VALUE][0:100]
        plt.subplot(num_rows, num_cols, i+1)
        plt.errorbar(time, curve, yerr=stderrs[conc_str][VALUE][0:100])

        # Define the parameters
        max_rfu = fitting.Parameter(np.max(curve))
        min_rfu = fitting.Parameter(np.min(curve))
        k = fitting.Parameter(1e-4)

        def exp_func(t):
            return (max_rfu() - min_rfu())*(1 - np.exp(-k()*t)) + min_rfu()
        result = fitting.fit(exp_func, [max_rfu, min_rfu, k], curve, time)
        cov_x = result[1]

        plt.plot(time, exp_func(time), color='r')
        plt.title(conc_str)
        plt.xlabel('Time (sec)')
        plt.ylabel('ANTS (RFU)')
        max_rfu_arr[i] = [max_rfu(), np.sqrt(cov_x[0, 0])]
        min_rfu_arr[i] = [min_rfu(), np.sqrt(cov_x[1, 1])]
        k_arr[i] = [k(), np.sqrt(cov_x[2, 2])]

    plt.tight_layout(pad=0, w_pad=0, h_pad=0)

    # Plot values of fitted parameters
    plt.figure()
    plt.errorbar(conc_list, k_arr[:,0], yerr=k_arr[:,1], marker='o')
    plt.figure()
    plt.errorbar(conc_list, min_rfu_arr[:,0], yerr=min_rfu_arr[:,1], marker='o')
    plt.figure()
    plt.errorbar(conc_list, max_rfu_arr[:,0], yerr=max_rfu_arr[:,1], marker='o')
    plt.figure()
    plt.plot(conc_list, (max_rfu_arr[:,0] - min_rfu_arr[:,0])/min_rfu_arr[:,0],
             marker='o')
    plt.title('Percent increase')
    import ipdb; ipdb.set_trace()

def approach_2(wells, layout):
    """In this approach we perform MCMC fitting of the replicate averages
    for each concentration and plot the uncertainty associated with each
    on the final titration plot."""
    # Get the average trajectories, with standard errors
    (avgs, stderrs) = averages(timecourse_wells, layout, stderr=True)

    iter=50000
    burn=iter/4.
    thin=50.
    num_samples = (iter - burn) / thin

    num_concs = len(layout.keys())
    min_rfu_arr = np.zeros((num_concs, 2))
    max_rfu_arr = np.zeros((num_concs, 2))
    k_arr = np.zeros((num_concs, 2))
    conc_list = np.zeros(num_concs)
    k_traces = np.zeros((num_concs, num_samples))

    def exp_func(max, min, k, time):
        return (max - min)*(1 - np.exp(-k * time)) + min

    plt.figure()
    num_rows = 4
    num_cols = 3
    for i, conc_str in enumerate(layout.keys()):
        print conc_str
        conc = float(conc_str.split(' ')[1])
        conc_list[i] = conc
        time = avgs[conc_str][TIME][0:100]
        curve = avgs[conc_str][VALUE][0:100]
        plt.subplot(num_rows, num_cols, i+1)
        plt.errorbar(time, curve, yerr=stderrs[conc_str][VALUE][0:100], color='r')

        max_rfu = pymc.Uniform('max_rfu', np.min(curve), 5*np.min(curve))
        min_rfu = pymc.Normal('min_rfu', mu=np.min(curve),
                              tau=1/((np.min(curve) * 0.05)**2))
        # Range between the min and max timescale
        # TODO Should get rid of magic numbers, calculate from data instead
        k = pymc.Uniform('k', 1 / 36000., 1 / 120.)

        @pymc.stochastic(plot=True, observed=True, verbose=0)
        def exp_fit(max=max_rfu, min=min_rfu, k=k, value=curve):
            f = exp_func(max, min, k, time)
            #cv = 0.01 # Magic number, should estimate from data
            #variance = (np.min(curve) * cv)**2
            variance = stderrs[conc_str][VALUE][:100] **2
            log_lkl = -np.sum(((curve - f)**2) / (2 * variance))
            return log_lkl

        mcmc = pymc.MCMC([max_rfu, min_rfu, k, exp_fit])
        mcmc.sample(iter=iter, burn=burn, thin=thin)
        print

        # TODO: should select randomly, not just last n samples
        num_to_plot = 500
        max_trace = mcmc.trace('max_rfu')[:]
        min_trace = mcmc.trace('min_rfu')[:]
        k_trace = mcmc.trace('k')[:]
        for j in range(num_to_plot):
            plt.plot(time, exp_func(max_trace[-j], min_trace[-j], k_trace[-j], time),
                     alpha=0.05, color='b')

        #pymc.Matplot.plot(mcmc)

        plt.title(conc_str)
        plt.xlabel('Time (sec)')
        plt.ylabel('ANTS (RFU)')
        #import ipdb; ipdb.set_trace()
        max_rfu_arr[i] = [np.mean(max_trace), np.std(max_trace)]
        min_rfu_arr[i] = [np.mean(min_trace), np.std(min_trace)]
        k_arr[i] = [np.mean(k_trace), np.std(k_trace)]
        k_traces[i, :] = k_trace

    plt.tight_layout(pad=0, w_pad=0, h_pad=0)

    # Plot values of fitted parameters
    plt.figure()
    plt.errorbar(conc_list, k_arr[:,0], yerr=k_arr[:,1])
    plt.figure()
    plt.errorbar(conc_list, min_rfu_arr[:,0], yerr=min_rfu_arr[:,1])
    plt.figure()
    plt.errorbar(conc_list, max_rfu_arr[:,0], yerr=max_rfu_arr[:,1])
    plt.figure()
    plt.boxplot(k_traces[:-1].T)

def approach_3(timecourse_wells, layout):
    """Global fit of all curves by least squares."""
    # Get the average trajectories, with standard errors
    (avgs, stderrs) = averages(timecourse_wells, layout, stderr=True)
    #(avgs, stderrs) = averages(timecourse_wells, layout, stderr=False)

    # Build full matrix of all timecourse averages and SDs
    num_concs = len(layout.keys())
    num_timepoints = 100
    time_arr = np.zeros((num_timepoints, num_concs))
    avg_arr = np.zeros((num_timepoints, num_concs))
    std_arr = np.zeros((num_timepoints, num_concs))
    conc_list = np.zeros(num_concs)
    for i, conc_str in enumerate(layout.keys()):
        conc = float(conc_str.split(' ')[1])
        conc_list[i] = conc
        time_arr[:, i] = avgs[conc_str][TIME][:num_timepoints]
        avg_arr[:, i] = avgs[conc_str][VALUE][:num_timepoints]
        std_arr[:, i] = stderrs[conc_str][VALUE][:num_timepoints]

    # Starting value for highest concentration
    min_initial = np.min(avg_arr[0,:])
    max_initial = avg_arr[0,0]
    max_conc = conc_list[0]
    min_final = np.min(avg_arr[-1,:])
    max_final = avg_arr[-1,0]
    min_int = fitting.Parameter(min_initial)
    min_slope = fitting.Parameter((max_initial - min_initial) / max_conc)
    max_int = fitting.Parameter(min_final)
    max_slope = fitting.Parameter((max_final - min_final) / max_conc)
    k = fitting.Parameter(1.2e-4)

    avg_arr_lin = avg_arr.T.reshape(num_timepoints * num_concs)
    def exp_func(t, conc, min_int, min_slope, max_int, max_slope, k):
        max_val = max_slope * conc + max_int
        min_val = min_slope * conc + min_int
        return (max_val - min_val) * (1 - np.exp(-k * t)) + min_val

    def exp_matrix(t):
        result = np.zeros((num_timepoints*num_concs))
        for i, conc in enumerate(conc_list):
            result[i*num_timepoints:(i+1)*num_timepoints] = \
                    exp_func(t, conc, min_int(), min_slope(),
                                   max_int(), max_slope(), k())
        return result

    fitting.fit(exp_matrix, [min_int, min_slope, max_int, max_slope, k],
                avg_arr_lin, time_arr[:,0])

    plt.figure()
    plt.plot(conc_list, avg_arr[0,:], marker='o')
    plt.plot(conc_list, [min_slope()*conc + min_int() for conc in conc_list], marker='o')
    plt.figure()
    num_rows = 4
    num_cols = 3
    def well_exp(initial_val, t):
        return 0.05*initial_val*(1 - np.exp(-1.5e-4 * t)) + initial_val

    for i, conc_str in enumerate(layout.keys()):
        print conc_str
        conc = float(conc_str.split(' ')[1])
        conc_list[i] = conc
        time = avgs[conc_str][TIME][0:100]
        curve = avgs[conc_str][VALUE][0:100]
        plt.subplot(num_rows, num_cols, i+1)
        plt.errorbar(time, curve, yerr=stderrs[conc_str][VALUE][0:100], color='b')
        plt.plot(time, exp_func(time, conc, min_int(), min_slope(),
                                max_int(), max_slope(), k()), color='r')
        plt.plot(time, well_exp(curve[0], time), color='g')
        plt.title(conc_str)

    import ipdb; ipdb.set_trace()


if __name__ == '__main__':
    plt.ion()

    #plot_data()
    #approach_3(timecourse_wells, layout)
    #approach_2(timecourse_wells, layout)
    #approach_1(timecourse_wells, layout)
    sys.exit()


    for well_name in layout['Liposomes 0.078 mg/ml']:
        curve = timecourse_wells[well_name][VALUE]
        time = timecourse_wells[well_name][TIME]
        plt.figure()
        plt.plot(time, curve)

        A = fitting.Parameter(450000)
        k = fitting.Parameter(1e-6)
        b = fitting.Parameter(390000)
        def exp_func(t): return A()*(1 - np.exp(-k()*t)) + b()
        fitting.fit(exp_func, [A, k, b], curve, time)
        plt.plot(time, exp_func(time))

    import pymc


    sys.exit()

    #import layout_131030 as ds # for dilution series
    #ds.num_concs = 12
    #ds.num_reps = 8
    #ds.num_timepoints = 139

    #(timecourse_wells, timecourse_averages, timecourse_stds,
    #        conc_list, means, stds, log_concs, log_means, log_stds) = \
    #                ds.get_wells('140221_ANTS_titration_timecourse.csv', layout,
    #                             reader='wallac')
    #ds.plot_bar_plot(means, stds)
    #ds.plot_cvs(means, stds, log_concs)
    #ds.plot_dilution_series(conc_list, means, stds, log_concs, log_means, log_stds)
    #ds.plot_fits(means, log_means, conc_list, log_concs)

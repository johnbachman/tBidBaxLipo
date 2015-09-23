import collections
import sys
import os
from copy import deepcopy
import pickle
import itertools

import numpy as np
from matplotlib import pyplot as plt
import emcee
from emcee.utils import MPIPool
import corner

from matplotlib.ticker import MultipleLocator, ScalarFormatter
from scipy.stats import linregress
from scipy.stats import distributions as dist

from tbidbaxlipo.util.plate_assay import *
import tbidbaxlipo.data
from tbidbaxlipo.util import fitting, set_fig_params_for_publication, emcee_fit
from tbidbaxlipo.plots import titration_fits as tf
from tbidbaxlipo.util import fitting, colors, set_fig_params_for_publication

#### LOAD THE DATA ####

# All wells with liposomes ~0.1 mg/mL
data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                              '141125_Bid_saturation_NBD_Bax_timecourse.txt'))

# Parse the data
timecourse_wells = read_flexstation_kinetics(timecourse_file)

# The raw layout will contain all wells in the dataset
raw_layout = collections.OrderedDict(
            reversed(dose_series_replicate_list('cBid', 800, 0.5, 12, True,
                                                'nM', lowest_first=False,
                                                start_row='A', end_row='C',
                                                start_col=1, end_col=12)))
raw_layout['Bim BH3'] = row_wells('D')
raw_layout['Bim BH3, No NBD-Bax'] = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6']
raw_layout['No NBD-Bax'] = ['E7', 'E8', 'E9', 'E10', 'E11', 'E12']

#### REMOVE OUTLIERS ####

# We remove a few outliers to make the "clean" list of wells
clean_layout = deepcopy(raw_layout)
# Remove E1, which has a break in the timecourse at ~7000 sec
clean_layout['Bim BH3, No NBD-Bax'].remove('E1')
# Remove B4, which has a break in the timecourse at ~1000 sec
clean_layout['cBid 100.0 nM'].remove('B4')

#### ACCOUNT FOR START TIME OFFSET ####

# Because the different rows were pipetted at different times, and the read
# didn't start right away but after a delay, we account for these delays here
# by offsetting the different wells by the appropriate delay.
def time_offset_vector():
    read_time = 200 # seconds (3:20)
    row_times = {'A': 0,
                 'B': 20,
                 'C': 40,
                 'D': 60,
                 'E': 80,
             }
    offset_dict = {}
    # New time vector is TIME + offset, where
    # offset = (read_time - row_time)
    for row in ['A', 'B', 'C', 'D', 'E']:
        for col in range(1, 13):
            well = '%s%s' % (row, col)
            offset_dict[well] = read_time - row_times[row]
    return offset_dict

def add_time_offsets(timecourse_wells, offset_vector):
    for well_name in timecourse_wells.keys():
        well = timecourse_wells[well_name]
        well[TIME] += offset_vector[well_name]

# Updates time vector of timecourse wells in place
add_time_offsets(timecourse_wells, time_offset_vector())

#### DO BACKGROUND SUBTRACTION ####

# Take the average of the Bim BH3 background wells (no NBD-Bax)
(timecourse_averages, timecourse_stds) = \
                                averages(timecourse_wells, clean_layout)

# Get the background vector for the Bim BH3 condition
bim_bh3_bg_arr = timecourse_averages['Bim BH3, No NBD-Bax'][VALUE]
# Subtract this background from the Bim BH3 wells
clean_bim_bh3_layout = extract(['Bim BH3'], clean_layout)
bim_bh3_well_names = itertools.chain(*clean_bim_bh3_layout.values())
bim_bh3_wells = extract(bim_bh3_well_names, timecourse_wells)
bgsub_bim_bh3_wells = subtract_background(bim_bh3_wells, bim_bh3_bg_arr)
(bgsub_bim_bh3_avgs, bgsub_bim_bh3_stds) = averages(bgsub_bim_bh3_wells,
                                                    clean_bim_bh3_layout)

# Get the background vector for the Bid conditions
bid_bg_arr = timecourse_averages['No NBD-Bax'][VALUE]
# Subtract this background from the Bid wells
# Get the well names with Bid, and extract these from the full set
bid_cond_names = [s for s in clean_layout.keys()
                  if s.split(' ')[0] == 'cBid']
clean_bid_layout = extract(bid_cond_names, clean_layout)
# Flatten the well list in the layout to get a list of all Bid well names
bid_well_names = itertools.chain(*clean_bid_layout.values())
# Get the well data itself for the Bid wells
bid_wells = extract(bid_well_names, timecourse_wells)
# Subtract the Bid background from the Bid wells
bgsub_bid_wells = subtract_background(bid_wells, bid_bg_arr)
(bgsub_bid_avgs, bgsub_bid_stds) = averages(bgsub_bid_wells, clean_bid_layout)

#### PLOTTING FUNCTIONS #### 

def plot_raw_data():
    """Plots the raw data and various transformations of it."""
    for condition in raw_layout.keys():
        figure()
        plot_all(extract(raw_layout[condition], timecourse_wells))
        title(condition)

def plot_clean_data():
    """Plots the clean data."""
    for condition in clean_layout.keys():
        figure()
        plot_all(extract(clean_layout[condition], timecourse_wells))
        title(condition)

def plot_averages():
    plt.figure()
    plot_all(bgsub_bid_avgs)
    plot_all(bgsub_bim_bh3_avgs)

def plot_bg():
    figure()
    tc_avg = timecourse_averages['Bim BH3, No NBD-Bax']
    tc_sd = timecourse_stds['Bim BH3, No NBD-Bax']
    errorbar(tc_avg[TIME], tc_avg[VALUE], label='Bim BH3 + Liposomes',
             yerr=tc_sd)

    tc_avg = timecourse_averages['No NBD-Bax']
    tc_sd = timecourse_stds['No NBD-Bax']
    errorbar(tc_avg[TIME], tc_avg[VALUE], label='Liposomes', yerr=tc_sd)
    legend(loc='lower right')

#### FITTING FUNCTIONS ####

nwalkers = 1000
burn_steps = 2000
sample_steps = 150
time = bgsub_bim_bh3_wells['D1'][TIME]
data = [bgsub_bim_bh3_wells['D%d' % i][VALUE]
        for i in range(1, 13)]
num_hyper_params = 6
ndim = 3 * len(data) + num_hyper_params

# Estimate SD of data by calculating SD from last 30 points of
# trajectory (I did this with full Bayesian estimation once and found
# it to be virtually identical to the value from the SD calculated from
# the points).
data_sd = np.std(data[0][-30:])
print "SD", data_sd

# Parameters will be in order (Fmax, k, F0)
def fit_func(position):
    (fmax, k, f0) = position
    return (fmax*f0) * (1 - np.exp(-k * time)) + f0

def likelihood(position):
    # Calculate the predicted timecourse
    loglkl = 0
    # The hyperparameters are at the end of the parameter array in the
    # order defined by the tuple below.
    hp_ix = len(data)*3
    (fmax_mean, fmax_sd, k_mean, k_sd, f0_mean, f0_sd) = position[hp_ix:]
    # Iterate over the trajectories in the dataset
    for d_ix, data_i in enumerate(data):
        # Get the parameters for this particular trajectory
        pos_i = position[(d_ix*3):(3*(d_ix+1))]
        # Get the predicted timecourse
        ypred = fit_func(pos_i)
        # This normalization factor is invariant so it shouldn't matter
        #err += -len(ypred) * np.log(data_sd * np.sqrt(2 * np.pi))
        # Calculate the log-likelihood
        loglkl += -np.sum((ypred - data_i) **2 / (2 * data_sd **2))
        (fmax, k, f0) = pos_i
        fmax_p = dist.norm.logpdf(fmax, loc=fmax_mean, scale=fmax_sd)
        k_p = dist.norm.logpdf(k, loc=k_mean, scale=k_sd)
        f0_p = dist.norm.logpdf(f0, loc=f0_mean, scale=f0_sd)
        if np.any(np.isnan([fmax_p, k_p, f0_p])):
            return -np.inf
        else:
            loglkl += fmax_p + k_p + f0_p
    return loglkl

def prior(position):
    """
    for d_ix in range(len(data)):
        (fmax, k, f0) = position[(d_ix*3):(3*(d_ix+1))]
        if fmax < 1 or fmax > 6:
            return -np.inf
        if k < 6e-5 or k > 1e-3:
            return -np.inf
        if f0 < 2 or f0 > 3:
            return -np.inf
    """
    hp_ix = len(data)*3
    (fmax_mean, fmax_sd, k_mean, k_sd, f0_mean, f0_sd) = position[hp_ix:]
    # Make sure all SDs are positive
    if fmax_sd < 0 or k_sd < 0 or f0_sd < 0:
        return -np.inf
    # The fmax should be between 1 and 6
    if fmax_mean < 1 or fmax_mean > 6:
        return -np.inf
    # k should be between 6e-5 and 1e-3
    if k_mean < 6e-5 or k_mean > 1e-3:
        return -np.inf
    # f0 should be between 2 and 3
    if f0_mean < 2 or f0_mean > 3:
        return -np.inf
    # Otherwise, the log-prior is 0 because all priors are uniform
    return 0.

def posterior(position):
    return likelihood(position) + prior(position)

def fit_bim_bh3_curves(p0=None):
    # Choose initial position
    if p0 is None:
        p0 = np.zeros((nwalkers, ndim))
        for walk_ix in range(nwalkers):
            for d_ix in range(len(data)):
                p0[walk_ix, d_ix*3] = np.random.uniform(1, 6)
                p0[walk_ix, d_ix*3 + 1] = np.random.uniform(6e-5, 1e-3)
                p0[walk_ix, d_ix*3 + 2] = np.random.uniform(2, 3)
            hp_ix = len(data)*3
            p0[walk_ix, hp_ix] = np.random.uniform(1,6) # fmax mean
            p0[walk_ix, hp_ix + 1] = np.random.uniform(0,1) # fmax sd
            p0[walk_ix, hp_ix + 2] = np.random.uniform(6e-5, 1e-3) # k mean
            p0[walk_ix, hp_ix + 3] = np.random.uniform(0,1e-1) # k sd
            p0[walk_ix, hp_ix + 4] = np.random.uniform(2,3) # f0 mean
            p0[walk_ix, hp_ix + 5] = np.random.uniform(0,1) # f0 sd

    #plt.figure()
    #for d_ix, data_i in enumerate(data):
    #    plt.plot(time, data_i, color=colors[d_ix])
    #    plt.plot(time, fit_func(p0[0, d_ix*3:(d_ix+1)*3]), color='k')

    # Initialize the MPI pool
    pool = MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    # Get the sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, posterior, pool=pool)
    # Burn-in
    print("Burn-in sampling...")
    pos, prob, state = sampler.run_mcmc(p0, burn_steps, storechain=False)
    sampler.reset() 
    # Main sampling
    print("Main sampling...")
    sampler.run_mcmc(pos, sample_steps)

    # Close the pool!
    pool.close()

    # Pickle the sampler
    sampler.pool = None
    with open('bimbh3_141125_2.pck','w') as f:
        pickle.dump(sampler, f)

    return sampler

def plot_chain(sampler):
    # Check convergence
    plt.figure()
    plt.plot(sampler.lnprobability.T)

    # Plot maximum likelihood
    plt.figure()
    ml_ix = np.unravel_index(np.argmax(sampler.lnprobability),
                             sampler.lnprobability.shape)
    ml_pos = sampler.chain[ml_ix]
    for d_ix, data_i in enumerate(data):
        plt.plot(time, data_i)
        plt.plot(time, fit_func(ml_pos[d_ix*3:(d_ix+1)*3]), color='k')

    # Plot sampling of trajectories parameters
    plt.figure()
    for d_ix, data_i in enumerate(data):
        plt.plot(time, data_i)
    num_plot_samples = 100
    num_tot_steps = sampler.flatchain.shape[0]
    for s_ix in range(num_plot_samples):
        p_ix = np.random.randint(num_tot_steps)
        p_samp = sampler.flatchain[p_ix]
        for d_ix in range(len(data)):
            plt.plot(time, fit_func(p_samp[d_ix*3:(d_ix+1)*3]), alpha=0.1,
                     color='k')

    # Triangle plots
    for d_ix in range(len(data)):
        corner.corner(sampler.flatchain[:,d_ix*3:(d_ix+1)*3])
    corner.corner(sampler.flatchain[:,3*len(data):])

def get_bim_bh3_curves():
    bim_bh3_curves = {}
    # Now iterate over each curve, and estimate the initial value by fitting
    # a line through the first 20 points and getting the intercept.
    figure()
    for well_name, well_value in bgsub_bim_bh3_wells.iteritems():
        t = well_value[TIME]
        v = well_value[VALUE]
        # Do the linear fit over the first numpts
        numpts = 20
        lin_fit = linregress(t[:numpts], v[:numpts])
        slope = lin_fit[0]
        intercept = lin_fit[1]
        intercept = 2.5
        norm_curve = v / float(intercept)
        # Build up the entry for the normalized curve
        bim_bh3_curves[well_name] = []
        bim_bh3_curves[well_name].append(t)
        bim_bh3_curves[well_name].append(norm_curve)
        # Plot the line fits over the data
        plot(t[:numpts], v[:numpts])
        line_pts = np.linspace(0, t[numpts-1], 20)
        plot(line_pts, line_pts * slope + intercept, 'k')

    return bim_bh3_curves

def fit_timecourses(layout, timecourses, plot=True):
    """Fit all timecourses in the experiment."""
    # We'll be doing a 3 parameter fit
    num_params = 3
    # This is where we'll store the results of the fitting
    param_dict = collections.OrderedDict()
    # Iterate over all of the conditions in the layout
    for cond_name, rep_list in layout.iteritems():
        # How many replicates for this condition?
        num_reps = len(rep_list)
        # Set up a plot to show the fits
        if plot:
            plt.figure(cond_name, figsize=(6, 2), dpi=150)
        plots_per_row = 3
        # Create an entry in the resulting param dict for this condition.
        # Results will be stored in a matrix with shape (num_reps, num_params)
        param_matrix = np.zeros((num_reps, num_params))
        # Now iterate over every well in the replicate list
        for rep_ix, well_name in enumerate(rep_list):
            # Get the time and fluorescence vectors
            time = timecourses[well_name][TIME]
            value = timecourses[well_name][VALUE]
            # Set up the fitting procedure: first define the parameters
            fmax = fitting.Parameter(2.1)
            f0 = fitting.Parameter(2.3)
            k = fitting.Parameter(1.4e-4)
            # Define the fitting function
            def exp_func(t):
                return (fmax() * f0()) * (1 - np.exp(-k() * t)) + f0()
            # Do the fit
            res = fitting.fit(exp_func, [fmax, k, f0], value, time)
            # Save the values in the order: (fmax, k, f0)
            param_matrix[rep_ix, 0] = fmax()
            param_matrix[rep_ix, 1] = k()
            param_matrix[rep_ix, 2] = f0()
            # Plot into the appropriate subplot
            if plot:
                plt.subplot(1, num_plots_per_row, rep_ix+1)
                plt.plot(time, value, color='r')
                plt.plot(time, exp_func(time), color='k', linewidth=2)
                plt.xlabel('Time (sec)')
                plt.ylabel('RFU')
        # Store the resulting parameter matrix in the dict entry for this cond
        param_dict[cond_name] = param_matrix

    return param_dict

def plot_k_titration(fits):
    bid_concs = []
    k_means = []
    k_ses = []

    for cond_name in bid_fits.keys():
        fits = bid_fits[cond_name]
        bid_conc = float(cond_name.split(' ')[1])
        bid_concs.append(bid_conc)
        num_reps = fits.shape[0]
        k_means.append(np.mean(fits[:,1]))
        k_ses.append(np.std(fits[:,1]) / np.sqrt(num_reps))

    bid_concs = np.array(bid_concs)
    k_means = np.array(k_means)

    plt.figure('k')
    plt.errorbar(np.log10(bid_concs), k_means, yerr=k_ses, color='k')

    kd1 = fitting.Parameter(1.)
    k0 = fitting.Parameter(2e-5)
    kmax1 = fitting.Parameter(2e-4)
    slope = fitting.Parameter(1e-8)

    def hill_func_line(concs):
        return (k0() +
               ((kmax1() * concs) / (kd1() + concs)) +
                (slope() * concs))
    res = fitting.fit(hill_func_line, [kd1, k0, kmax1, slope],
                      np.array(k_means), np.array(bid_concs))
    plt.figure('k')
    plt.plot(np.log10(bid_concs), hill_func_line(bid_concs), color='b')
    #plt.plot(np.log10(bid_concs), slope() * bid_concs, color='r')
    #plt.plot(np.log10(bid_concs),
    #         k0() + ((kmax1()*bid_concs) / (kd1() + bid_concs)), color='g')

    kd1 = fitting.Parameter(1.)
    k0 = fitting.Parameter(7e-5)
    kmax1 = fitting.Parameter(2e-4)

    def hill_func_k(concs):
        return (k0() + ((kmax1() * concs) / (kd1() + concs)))

    res = fitting.fit(hill_func_k, [kd1, k0, kmax1],
                      np.array(k_means), np.array(bid_concs))
    plt.figure('k')
    plt.plot(np.log10(bid_concs), hill_func_k(bid_concs), color='r')


if __name__ == '__main__':
    plt.ion()
    # plot_raw_data()
    # plot_clean_data()
    # plot_bg()
    # plot_averages()
    # get_bim_bh3_curves()
    #with open('bimbh3_141125.pck') as f:
    #    sampler = pickle.load(f)
    #sampler = fit_bim_bh3_curves(p0=sampler.chain[:,-1,:])
    #plot_chain(sampler)
    bid_fits = fit_timecourses(clean_bid_layout, bgsub_bid_wells, plot=False)
    plot_k_titration(bid_fits)

    """

    bid_concs = []
    fmax_means = []
    fmax_ses = []
    f0_means = []
    f0_ses = []

    for cond_name in bid_fits.keys():
        fits = bid_fits[cond_name]
        bid_conc = float(cond_name.split(' ')[1])
        bid_concs.append(bid_conc)
        num_reps = fits.shape[0]
        fmax_means.append(np.mean(fits[:,0]))
        fmax_ses.append(np.std(fits[:,0]) / np.sqrt(num_reps))
        f0_means.append(np.mean(fits[:,2]))
        f0_ses.append(np.std(fits[:,2]) / np.sqrt(num_reps))

    bid_concs = np.array(bid_concs)
    fmax_means = np.array(fmax_means)

    plt.figure('fmax')
    plt.errorbar((bid_concs), fmax_means, yerr=fmax_ses, linestyle='')
    plt.figure('f0')
    plt.errorbar((bid_concs), f0_means, yerr=f0_ses, linestyle='')

    fmaxd = fitting.Parameter(1.)
    fmax0 = fitting.Parameter(0.7)
    fmaxmax = fitting.Parameter(2.2)

    def hill_func_fmax(concs):
        return fmax0() + ((fmaxmax() * concs) / (fmaxd() + concs))

    res = fitting.fit(hill_func_fmax, [fmaxd, fmax0, fmaxmax],
                      np.array(fmax_means), np.array(bid_concs))
    plt.figure('fmax')
    plt.plot(bid_concs, hill_func_fmax(bid_concs), color='k')
    """

import collections
import sys
import os
from copy import deepcopy

import numpy as np
from matplotlib import pyplot as plt
import emcee
import triangle

from matplotlib.ticker import MultipleLocator, ScalarFormatter
from scipy.stats import linregress
from scipy.stats import distributions as dist

from tbidbaxlipo.util.plate_assay import *
import tbidbaxlipo.data
from tbidbaxlipo.util import fitting, set_fig_params_for_publication, emcee_fit
from tbidbaxlipo.plots import titration_fits as tf
from tbidbaxlipo.util import fitting, colors, set_fig_params_for_publication

# All wells with liposomes ~0.1 mg/mL
# Load the data
data_path = os.path.dirname(sys.modules['tbidbaxlipo.data'].__file__)
timecourse_file = os.path.abspath(os.path.join(data_path,
                              '141125_Bid_saturation_NBD_Bax_timecourse.txt'))

# Parse the data
timecourse_wells = read_flexstation_kinetics(timecourse_file)

# The raw layout will contain all wells in the dataset
raw_layout = collections.OrderedDict(
            dose_series_replicate_list('cBid', 800, 0.5, 12, True, 'nM',
                                    lowest_first=False,
                                    start_row='A', end_row='C',
                                    start_col=1, end_col=12))
raw_layout['Bim BH3'] = row_wells('D')
raw_layout['Bim BH3, No NBD-Bax'] = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6']
raw_layout['No NBD-Bax'] = ['E7', 'E8', 'E9', 'E10', 'E11', 'E12']

# We remove a few outliers to make the "clean" list of wells
clean_layout = deepcopy(raw_layout)
# Remove E1, which has a break in the timecourse at ~7000 sec
clean_layout['Bim BH3, No NBD-Bax'].remove('E1')
# Remove B4, which has a break in the timecourse at ~1000 sec
clean_layout['cBid 100.0 nM'].remove('B4')

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

# Take the average of the Bim BH3 background wells (no NBD-Bax)
(timecourse_averages, timecourse_stds) = \
                                averages(timecourse_wells, clean_layout)

# Get the background vector for the Bim BH3 condition
bim_bh3_bg_arr = timecourse_averages['Bim BH3, No NBD-Bax'][VALUE]
# Subtract this background from the Bim BH3 wells
bim_bh3_wells = extract(clean_layout['Bim BH3'], timecourse_wells)
bgsub_bim_bh3_wells = subtract_background(bim_bh3_wells, bim_bh3_bg_arr)

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

def fit_bim_bh3_curves():
    nwalkers = 1000
    burn_steps = 1000
    sample_steps = 50
    tc1 = bgsub_bim_bh3_wells['D1']
    time = tc1[TIME]
    data1 = tc1[VALUE]
    tc2 = bgsub_bim_bh3_wells['D2']
    data2 = tc2[VALUE]
    tc3 = bgsub_bim_bh3_wells['D3']
    data3 = tc3[VALUE]
    tc4 = bgsub_bim_bh3_wells['D4']
    data4 = tc4[VALUE]
    tc5 = bgsub_bim_bh3_wells['D5']
    data5 = tc5[VALUE]
    tc6 = bgsub_bim_bh3_wells['D6']
    data6 = tc6[VALUE]
    tc7 = bgsub_bim_bh3_wells['D7']
    data7 = tc7[VALUE]
    tc8 = bgsub_bim_bh3_wells['D8']
    data8 = tc8[VALUE]
    tc9 = bgsub_bim_bh3_wells['D9']
    data9 = tc9[VALUE]
    tc10 = bgsub_bim_bh3_wells['D10']
    data10 = tc10[VALUE]
    tc11 = bgsub_bim_bh3_wells['D11']
    data11 = tc11[VALUE]
    tc12 = bgsub_bim_bh3_wells['D12']
    data12 = tc12[VALUE]
    data = [data1, data2, data3, data4, data5, data6,
            data7, data8, data9, data10, data11, data12]
    num_hyper_params = 2
    ndim = 3 * len(data) + num_hyper_params

    # Estimate SD by calculating SD from last 30 points of trajectory
    data_sd = np.std(data1[-30:])
    print "SD", data_sd

    # Parameters will be in order (Fmax, k, F0)
    def fit_func(position):
        (fmax, k, f0) = position
        return (fmax*f0) * (1 - np.exp(-k * time)) + f0

    def likelihood(position):
        # Calculate the predicted timecourse
        err = 0
        fmax_ix = len(data)*3
        (fmax_mean, fmax_sd) = position[fmax_ix:fmax_ix+2]
        for d_ix, data_i in enumerate(data):
            pos_i = position[(d_ix*3):(3*(d_ix+1))]
            ypred = fit_func(pos_i)
            err += -len(ypred) * np.log(data_sd * np.sqrt(2 * np.pi))
            err += -np.sum((ypred - data_i) **2 / (2 * data_sd **2))
            (fmax, k, f0) = pos_i
            fmax_p = dist.norm.logpdf(fmax, loc=fmax_mean, scale=fmax_sd)
            if np.isnan(fmax_p):
                return -np.inf
            else:
                err += fmax_p
        return err

    def prior(position):
        for d_ix in range(len(data)):
            (fmax, k, f0) = position[(d_ix*3):(3*(d_ix+1))]
            if fmax < 1 or fmax > 6:
                return -np.inf
            if k < 6e-5 or k > 1e-3:
                return -np.inf
            if f0 < 2 or f0 > 3:
                return -np.inf
        fmax_ix = len(data)*3
        (fmax_mean, fmax_sd) = position[fmax_ix:fmax_ix+2]
        if fmax_sd < 0:
            return -np.inf
        if fmax_mean < 1 or fmax_mean > 6:
            return -np.inf
        return 0.

    def posterior(position):
        return likelihood(position) + prior(position)

    # Choose initial position
    p0 = np.zeros((nwalkers, ndim))
    for walk_ix in range(nwalkers):
        for d_ix in range(len(data)):
            p0[walk_ix, d_ix*3] = np.random.uniform(1, 6)
            p0[walk_ix, d_ix*3 + 1] = np.random.uniform(6e-5, 1e-3)
            p0[walk_ix, d_ix*3 + 2] = np.random.uniform(2, 3)
        fmax_ix = len(data)*3
        p0[walk_ix, fmax_ix] = np.random.uniform(1,6)
        p0[walk_ix, fmax_ix + 1] = np.random.uniform(0,1)

    plt.figure()
    for d_ix, data_i in enumerate(data):
        plt.plot(time, data_i, color=colors[d_ix])
        plt.plot(time, fit_func(p0[0, d_ix*3:(d_ix+1)*3]), color='k')

    # Get the sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, posterior)
    # Burn-in
    print("Burn-in sampling...")
    pos, prob, state = sampler.run_mcmc(p0, burn_steps, storechain=False)
    sampler.reset() 
    # Main sampling
    print("Main sampling...")
    sampler.run_mcmc(pos, sample_steps)

    # Check convergence
    plt.figure()
    plt.plot(sampler.lnprobability.T)

    # Plot maximum likelihood
    plt.figure()
    ml_ix = np.unravel_index(np.argmax(sampler.lnprobability),
                             sampler.lnprobability.shape)
    ml_pos = sampler.chain[ml_ix]
    for d_ix, data_i in enumerate(data):
        plt.plot(time, data_i, color=colors[d_ix])
        plt.plot(time, fit_func(ml_pos[d_ix*3:(d_ix+1)*3]), color='k')

    # Plot sampling of parameters
    plt.figure()
    for d_ix, data_i in enumerate(data):
        plt.plot(time, data_i, color=colors[d_ix])
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
        triangle.corner(sampler.flatchain[:,d_ix*3:(d_ix+1)*3])
    triangle.corner(sampler.flatchain[:,3*len(data):])

    import ipdb; ipdb.set_trace()
    np.mean(sampler.flatchain[:,:2:])
    return sampler

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

if __name__ == '__main__':
    plt.ion()
    sampler = fit_bim_bh3_curves()


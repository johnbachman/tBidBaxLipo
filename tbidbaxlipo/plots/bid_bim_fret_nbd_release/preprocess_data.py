import sys
import numpy as np
import pandas as pd
from copy import deepcopy
from tbidbaxlipo.data.parse_bid_bim_fret_nbd_release import df, nbd_residues
from tbidbaxlipo.data.parse_bid_bim_nbd_release import labeling_ratios
from tbidbaxlipo.util.calculate_error_variance import calc_err_var
from tbidbaxlipo.util.find_outliers import find_outliers
this_module = sys.modules[__name__]

# There are a number of outliers in the FRET data that Justin attributes
# to debris floating in the cuvette that produces momentary spikes. Here
# we explicitly remove several of these from the 54C data, setting them to
# NaN.
df_pre = deepcopy(df)

def remove_outliers(key, outliers):
    for outlier_ix in outliers:
        df_pre[key][outlier_ix] = np.nan

def remove_all_fret_outliers():
    # Bid/3/1: though negative, no terrible outliers
    # Bid/3/2: no terrible outliers
    # Bid 3/3: no terrible outliers
    # Bim/3/1
    remove_outliers(('Bim', 'FRET', '3', 1, 'VALUE'), [55, 58, 94, 105])
    # Bim/3/2
    remove_outliers(('Bim', 'FRET', '3', 2, 'VALUE'),
                    [12, 18, 20, 21, 30, 48, 56, 61, 69, 100, 106, 111, 114])
    # Bim/3/3: not too bad
    # Bid/15/1: goes negative, but not bad
    # Bid/15/2: not too bad
    # Bid/15/3: not too bad
    # Bim/15/1
    remove_outliers(('Bim', 'FRET', '15', 1, 'VALUE'), [39, 41, 43, 61, 76, 84])
    # Bim/15/2
    remove_outliers(('Bim', 'FRET', '15', 2, 'VALUE'), [3, 7, 8, 11, 13])
    # Bim/15/3
    remove_outliers(('Bim', 'FRET', '15', 3, 'VALUE'), [84])
    # Bid/36/1: pretty good
    # Bid/36/2: pretty good
    # Bid/36/3: pretty good
    # Bim/36/1
    remove_outliers(('Bim', 'FRET', '36', 1, 'VALUE'),
                    [21, 44, 59, 61, 63, 75, 82, 95, 96])
    # Bim/36/2
    # Bim/36/3
    remove_outliers(('Bim', 'FRET', '36', 3, 'VALUE'), [0, 28, 68, 97, 100])
    # Bid/47/1: two early outliers
    remove_outliers(('Bid', 'FRET', '47', 1, 'VALUE'), [7, 15])
    # Bid/47/2: Several outlying negative pts:
    remove_outliers(('Bid', 'FRET', '47', 2, 'VALUE'), [86, 97, 110, 111, 116])
    # Bid/47/3: No outliers
    # Bim/47/1: no outliers
    # Bim/47/2
    remove_outliers(('Bim', 'FRET', '47', 2, 'VALUE'), 
                    [1, 13, 20, 22, 24, 25, 54, 55, 57, 65, 68, 70, 71, 72, 76,
                     97, 117])
    # Bim/47/3
    remove_outliers(('Bim', 'FRET', '47', 3, 'VALUE'), [33, 39, 84])
    # Bid/54/1
    # Bid/54/2
    remove_outliers(('Bid', 'FRET', '54', 2, 'VALUE'), [111])
    # Bid/54/3
    # Bim/54/1: no outliers
    # Bim/54/2: no outliers
    # Bim/54/3: two bad negative ones, plus several more dubious ones
    remove_outliers(('Bim', 'FRET', '54', 3, 'VALUE'), [10, 21])
    # Bid/62/1: no outliers
    # Bid/62/2: many outliers, some worse than others
    remove_outliers(('Bid', 'FRET', '62', 2, 'VALUE'),
                    [23, 32, 44, 53, 56, 68, 72, 81, 83, 98, 112, 113])
    # Bid/62/3: this one is in rough shape. Nearly the entire trajectory is
    # negative, with only major outlying spikes above 0. This removes the spikes
    # but doesn't deal with the negative problem.
    remove_outliers(('Bid', 'FRET', '62', 3, 'VALUE'), 
                    [0, 1, 3, 11, 12, 15, 20, 29, 31, 33, 36, 37, 41, 44, 51,
                     53, 61, 64, 85, 86, 94, 96, 97, 98, 100, 111])
    # Bim/62/1
    remove_outliers(('Bim', 'FRET', '62', 1, 'VALUE'), [110])
    # Bim/62/2: this one is almost unusable, with negative values as low as
    # -169% FRET! Since most of the consistent points are positive, the simplest
    # heuristic is to just throw away the negative values.
    remove_outliers(('Bim', 'FRET', '62', 2, 'VALUE'),
                    [1, 5, 7, 13, 17, 24, 25, 26, 27, 29, 35, 36, 38,
                    42, 43, 46, 47, 48, 50, 51, 53, 55, 60, 65, 66, 69,
                    70, 77, 78, 82, 88, 90, 94, 95, 97, 100, 104, 105, 106,
                    107, 108, 109, 113, 115, 118])
    # Bim/62/3
    remove_outliers(('Bim', 'FRET', '62', 3, 'VALUE'),
                    [30, 62, 85, 91, 92, 101])
    # Bid/122/1
    remove_outliers(('Bid', 'FRET', '122', 1, 'VALUE'), [9])
    # Bid/122/2: pretty good
    # Bid/122/3: pretty good
    # Bim/122/1: pretty good
    # Bim/122/2
    remove_outliers(('Bim', 'FRET', '122', 2, 'VALUE'),
                    [7, 8, 28, 45, 84, 89, 90, 95, 96])
    # Bim/122/3
    remove_outliers(('Bim', 'FRET', '122', 3, 'VALUE'), [70, 71, 84, 86, 96])
    # Bid/126/1: OK
    # Bid/126/2: OK
    # Bid/126/3: OK
    # Bim/126/1
    remove_outliers(('Bim', 'FRET', '126', 1, 'VALUE'), [49])
    # Bim/126/2: mostly negative, with few positive outliers
    remove_outliers(('Bim', 'FRET', '126', 2, 'VALUE'), 
                    [2, 3, 7, 13, 15, 25, 29, 31, 37, 40, 43, 47, 57, 71, 74,
                    76, 78, 82, 102])
    # Bim/126/3
    remove_outliers(('Bim', 'FRET', '126', 3, 'VALUE'), [48, 49])
    # Bid/138/1
    remove_outliers(('Bid', 'FRET', '138', 1, 'VALUE'),
                    [7, 13, 14, 36, 47, 64, 81, 99, 100, 106, 107, 110, 113])
    # Bid/138/2: this one in rough shape, many drops in signal, hard
    # to draw a line for what to take out
    remove_outliers(('Bid', 'FRET', '138', 2, 'VALUE'),
                    [2, 19, 24, 40, 43, 47, 52, 74, 80, 81, 82, 85, 86, 89, ])
    # Bid/138/3
    remove_outliers(('Bid', 'FRET', '138', 3, 'VALUE'),
                    [7, 8, 9, 13, 31, 42, 51, 54, 55, 58, 73])
    # Bim/138/1
    remove_outliers(('Bim', 'FRET', '138', 1, 'VALUE'),
                    [53, 56, 62, 75, 90, 103, 114])
    # Bim/138/2
    remove_outliers(('Bim', 'FRET', '138', 2, 'VALUE'), [41])
    # Bim/138/3
    remove_outliers(('Bim', 'FRET', '138', 3, 'VALUE'), [32])
    # Bid/151/1
    remove_outliers(('Bid', 'FRET', '151', 1, 'VALUE'),
                    [2, 6, 76, 83, 88, 99, 103, 106, 118])
    # Bid/151/2: noise, but no outliers
    # Bid/151/3: noise, but no outliers
    # Bim/151/1
    remove_outliers(('Bim', 'FRET', '151', 1, 'VALUE'), [10])
    # Bim/151/2
    remove_outliers(('Bim', 'FRET', '151', 2, 'VALUE'),
                    [1, 50, 64, 87, 89, 94, 101, 103])
    # Bim/151/3: noisy, but no major outliers
    # Bid/175/1: very noisy, but hard to pick out outliers
    remove_outliers(('Bid', 'FRET', '175', 1, 'VALUE'),
                    [32, 41, 42, 51, 93, 107])
    # Bid/175/2: pretty good
    # Bid/175/3: pretty good
    # Bim/175/1: not too bad
    # Bim/175/2
    remove_outliers(('Bim', 'FRET', '175', 2, 'VALUE'),
                    [55, 57, 68, 71, 75, 87, 88, 106, 108, 114])
    # Bim/175/3
    remove_outliers(('Bim', 'FRET', '175', 3, 'VALUE'),
                    [44, 56, 57, 67, 69, 70, 72, 78])
    # Bid/184/1: pretty good
    # Bid/184/2
    remove_outliers(('Bid', 'FRET', '184', 2, 'VALUE'),
                    [26, 49, 52, 85, 94])
    # Bid/184/3
    remove_outliers(('Bid', 'FRET', '184', 3, 'VALUE'),
                    [2, 4, 39])
    # Bim/184/1: pretty good
    # Bim/184/2: pretty good
    # Bim/184/3: pretty good

activators = ['Bid', 'Bim']
observables = ['NBD', 'FRET']

#remove_all_fret_outliers()
for activator in activators:
    for nbd_residue in nbd_residues:
        if nbd_residue == 'WT':
            continue
        for rep_ix, rep_num in enumerate(range(1, 4)):
            key = (activator, 'FRET', nbd_residue, rep_num)
            time = df[key + ('TIME',)].values
            values = df[key + ('VALUE',)].values
            filt = find_outliers(values)
            df_pre[key + ('VALUE',)] = filt

# We're going to store the highest observed NBD fluorescence value as a
# reference to evaluate plausibility of fitted fluorescence values
max_nbd_value = -np.inf

# Instead of creating a list of ~60 variables for every replicate of every
# residue, we add entries to the __dict__ variable for this namespace,
# making this module appear to contain all ~60 variables. This way the data
# can be simply accessed programmatically as a variable, rather than having
# to modify the run_pt script to call a function.
for activator in activators:
    for nbd_residue in nbd_residues:
        # Skip the WT residue, since there's no NBD timecourse for it
        if nbd_residue == 'WT':
            continue
        # Get the number of timepoints based on first obs/rep;
        # If this turns out to not match, we will error later
        num_tpts = len(df[(activator, observables[0],
                           nbd_residue, 1, 'TIME')].values)
        # To estimate the experimental error, we fit each timecourse to a
        # polynomial equation and return the residuals; after processing all
        # replicates, we pool the residuals and calculate the standard error
        num_residuals = 50
        data_residuals = np.zeros((1, len(observables), num_residuals))
        #residuals_reps = []
        # Enumerate over the replicates
        for rep_ix, rep_num in enumerate(range(1, 4)):
            # We use the same timepoints for all observables
            time_var = df[(activator, observables[0],
                           nbd_residue, rep_num, 'TIME')].values
            # The fitting procedure expects a three-dimensional array for
            # each entry; the first dimension is reserved for multiple
            # concentrations, as in a titration, which we don't have; the
            # second is reserved for multiple observables, and the third
            # contains the timepoints.
            data_var = np.zeros((1, len(observables), len(time_var)))
            # Iterate over the observables
            for obs_ix, obs in enumerate(observables):
                timecourse = \
                     df[(activator, obs, nbd_residue, rep_num, 'VALUE')].values
                data_var[0, obs_ix, :] = timecourse

                # Get the residuals for this observable/rep
                (residuals, fig) = \
                   calc_err_var(timecourse, last_n_pts=num_residuals,
                                fit_type='cubic', plot=False)
                # Add the residuals for this rep to the list
                data_residuals[0, obs_ix, :] = residuals

                # NBD-specific variables
                if obs == 'NBD':
                    # Initial fluorescence value
                    f0 = timecourse[0]
                    # Update maximum fluorescence value
                    max_nbd_for_rep = np.max(timecourse)
                    if max_nbd_for_rep > max_nbd_value:
                        max_nbd_value = max_nbd_for_rep
                    # Get the minimum fluorescence allowable for this rep:
                    # 0.1 * the initial value
                    lbound = f0 * 0.1
                    # Get the maximum allowable fluorescence for this rep: the
                    # lesser of 10 times the initial value or a corrected
                    # absolute NBD fluorescence of 300k
                    max_nbd_fluorescence = 300000
                    mutant_absolute_max = \
                            labeling_ratios[nbd_residue] * max_nbd_fluorescence
                    mutant_relative_max = 10 * f0
                    ubound = min(mutant_absolute_max, mutant_relative_max)
                    # Upper and lower bounds for fluorescence parameters
                    data_nbd_var_name = 'data_%s_%s_%s_r%s' % \
                                (activator, obs, nbd_residue, rep_num)
                    lbound_name = data_nbd_var_name + '_lbound'
                    ubound_name = data_nbd_var_name + '_ubound'
                    f0_name = data_nbd_var_name + '_f0'
                    setattr(this_module, lbound_name, lbound)
                    setattr(this_module, ubound_name, ubound)
                    setattr(this_module, f0_name, f0)
                # -- end NBD specific block
            # -- end observables iteration
            # Variable names
            time_var_name = 'time_%s_%s_r%s' % \
                            (activator, nbd_residue, rep_num)
            data_var_name = 'data_%s_%s_r%s' % \
                            (activator, nbd_residue, rep_num)
            # Add these as module-level variables
            setattr(this_module, time_var_name, time_var)
            setattr(this_module, data_var_name, data_var)
        # -- end replicates iteration
        # Now that we've got the residuals for all reps, pool them and
        # calculate the standard deviation
        #pooled_residuals = np.array(residuals_reps).flatten()
        #pooled_err = np.std(pooled_residuals, ddof=1)
        # For the error estimates, the fitting procedure expects a 2D array,
        # with potentially distinct error estimates for each condition
        # (first axis) and observable (second axis)
        data_sigma_var = np.std(data_residuals, axis=2)
        # Add a module-level variable for the experimental error
        data_sigma_var_name = 'data_sigma_%s_%s' % (activator, nbd_residue)
        setattr(this_module, data_sigma_var_name, data_sigma_var)
    # -- end residues iteration
# -- end activators iteration


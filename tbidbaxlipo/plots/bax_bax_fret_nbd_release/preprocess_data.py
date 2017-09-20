import sys
import numpy as np
from copy import deepcopy
from tbidbaxlipo.data.parse_bax_bax_fret_nbd_release import df, nbd_residues
#from tbidbaxlipo.data.parse_bid_bim_nbd_release import labeling_ratios
from tbidbaxlipo.util.calculate_error_variance import calc_err_var
from tbidbaxlipo.util.find_outliers import find_outliers

this_module = sys.modules[__name__]

# There are a number of outliers in the FRET data that Justin attributes
# to debris floating in the cuvette that produces momentary spikes. Here
# we explicitly remove these, setting them to NaN.
df_pre = deepcopy(df)

activators = ['Bid']
observables = ['NBD', 'FRET']

# Remove FRET outliers
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
                    #mutant_absolute_max = \
                    #        labeling_ratios[nbd_residue] * max_nbd_fluorescence
                    # FIXME Get labeling ratios for these mutants?
                    mutant_absolute_max = max_nbd_fluorescence
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


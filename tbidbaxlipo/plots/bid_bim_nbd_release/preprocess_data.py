import sys
import numpy as np
from tbidbaxlipo.data.parse_bid_bim_nbd_release import df, nbd_residues
from copy import deepcopy

this_module = sys.modules[__name__]

# Copy the dataset before removing outliers
#df_pre = deepcopy(df)

# Remove outliers
pass

# Instead of creating a list of ~60 variables for every replicate of every
# residue, we add entries to the __dict__ variable for this namespace,
# making this module appear to contain all ~60 variables. This way the data
# can be simply accessed programmatically as a variable, rather than having
# to modify the run_pt script to call a function.
activator = 'Bid'
obs = 'NBD'
for nbd_residue in nbd_residues:
    # Skip the WT residue, since there's no NBD timecourse for it
    if nbd_residue == 'WT':
        continue
    # Enumerate over the replicates
    for rep_ix, rep_num in enumerate(range(1, 4)):
        time_var = df[(activator, obs, nbd_residue, rep_num, 'TIME')].values
        # The fitting procedure expects a three-dimensional array for each
        # entry; the first dimension is reserved for multiple concentrations,
        # as in a titration, which we don't have; the second is reserved for
        # multiple observables, which we also don't have for multiconf
        # fitting, and the third contains the timepoints.
        data_var = np.zeros((1, 1, len(time_var)))
        data_var[0, 0, :] = \
                df[(activator, obs, nbd_residue, rep_num, 'VALUE')].values
        # For the error estimates, the fitting procedure expects a 2-D array,
        # with potentially distinct error estimates for each condition
        # (first axis) and observable (second axis)
        data_sigma_var = np.array([[0.021085]])
        # Variable names
        time_var_name = 'time_%s_%s_%s_r%s' % \
                        (activator, obs, nbd_residue, rep_num)
        data_var_name = 'data_%s_%s_%s_r%s' % \
                        (activator, obs, nbd_residue, rep_num)
        data_sigma_var_name = 'data_sigma_%s_%s_%s_r%s' % \
                        (activator, obs, nbd_residue, rep_num)
        # Add these as module-level variables
        setattr(this_module, time_var_name, time_var)
        setattr(this_module, data_var_name, data_var)
        setattr(this_module, data_sigma_var_name, data_sigma_var)


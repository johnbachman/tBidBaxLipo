"""
A script for submitting multiple MCMC jobs to LSF.  Allows fitting multiple
models to multiple combinations of data trajectories with multiple combinations
of model observables.
"""

from tbidbaxlipo.mcmc import nbd_mcmc
import subprocess
import numpy as np

# Use this if you want to try fitting every model
# model_names = nbd_mcmc.model_names

# Otherwise, enumerate the models to fit here
model_names = ['tardt']
"""The models to fit to the data."""

nsteps = 100000
"""The number of steps in each chain."""

num_chains = 10
"""The number of chains to run for each model."""

queue = 'short'
"""The LSF queue to submit the jobs to."""

nbd_site_combos = [['c62']]
"""The combinations of NBD sites from the data to attempt to fit
simultaneously.  If an entry is a tuple with a single element (e.g.,
``('c3')``), fits only that particular trajectory.  """

nbd_obs_combos = [['iBax']]
"""The combinations of observables in the model to attempt to fit to the
corresponding combinations of NBD sites in the list ``nbd_site_combos``."""

T_inits = 10 ** np.linspace(0,7,15)
"""The initial temperatures at which to run the chains."""

thermo_temps = np.linspace(-7, 0, 15)
"""The initial temperatures at which to run the chains."""

if __name__ == '__main__':
    cmd_list = []
    for model in model_names:
        for nbd_site_combo in nbd_site_combos:
            nbd_site_str = '-'.join(nbd_site_combo)
            for nbd_obs_combo in nbd_obs_combos:
                nbd_obs_str = '-'.join(nbd_obs_combo)
                #for T_init in T_inits:
                for thermo_temp in thermo_temps:
                    for i in range(0, num_chains):
                        output_filename = '%s_%s_%s_%d_thermo%f_s%d.out' % \
                            (model, nbd_site_str, nbd_obs_str, nsteps, 
                                    thermo_temp, i)

                        cmd_list = ['bsub', '-W', '12:00',
                                    '-q', queue,
                                    '-o', output_filename,
                                    'python', '-m', 'tbidbaxlipo.mcmc.nbd_mcmc_run',
                                    'random_seed=%d' % i,
                                    "model=%s" % model,
                                    "nbd_sites=%s" % nbd_site_str,
                                    "nbd_observables=%s" % nbd_obs_str,
                                    "nsteps=%d" % nsteps,
                                    "thermo_temp=%.1f" % thermo_temp]
                        print ' '.join(cmd_list)
                        subprocess.call(cmd_list)


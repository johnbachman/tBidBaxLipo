"""
A script for submitting multiple MCMC jobs to LSF.  Allows fitting multiple
models to multiple combinations of data trajectories with multiple combinations
of model observables.
"""

from tbidbaxlipo.mcmc import nbd_mcmc
import subprocess

# Use this if you want to try fitting every model
# model_names = nbd_mcmc.model_names

# Otherwise, enumerate the models to fit here
model_names = ['tairdt']
"""The models to fit to the data."""

nsteps = 150000
"""The number of steps in each chain."""

num_chains = 10
"""The number of chains to run for each model."""

queue = 'short'
"""The LSF queue to submit the jobs to."""

nbd_site_combos = [('c3', 'c62')]
"""The combinations of NBD sites from the data to attempt to fit
simultaneously.  If an entry is a tuple with a single element (e.g.,
``('c3')``), fits only that particular trajectory.  """

nbd_obs_combos = [('iBax', 'Baxbh3'), ('iBax', 'Bax2')]
"""The combinations of observables in the model to attempt to fit to the
corresponding combinations of NBD sites in the list ``nbd_site_combos``."""

if __name__ == '__main__':
    cmd_list = []
    for model in model_names:
        for nbd_site_combo in nbd_site_combos:
            nbd_site_str = '-'.join(nbd_site_combo)
            for nbd_obs_combo in nbd_obs_combos:
                nbd_obs_str = '-'.join(nbd_obs_combo)
                for i in range(0, num_chains):
                    output_filename = '%s_%s_%s_%d_s%d.out' % \
                        (model, nbd_site_str, nbd_obs_str, nsteps, i)

                    cmd_list = ['bsub', '-W', '12:00',
                                '-q', queue,
                                '-o', output_filename,
                                'python', '-m', 'tbidbaxlipo.mcmc.nbd_mcmc_run',
                                'random_seed=%d' % i,
                                "model=%s" % model,
                                "nbd_sites=%s" % nbd_site_str,
                                "nbd_observables=%s" % nbd_obs_str,
                                "nsteps=%d" % nsteps]
                    print ' '.join(cmd_list)
                    #subprocess.call(cmd_list)


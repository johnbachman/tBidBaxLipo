"""
A script for submitting multiple MCMC jobs for the NBD plate data to LSF using
MPI-based parallel tempering.  Allows fitting models with alternative numbers
of assumed conformational states to different replicates for each NBD-labeled
Bax mutant.
"""

from tbidbaxlipo.mcmc import nbd_plate_mcmc

import subprocess
import numpy as np

# Use this if you want to try fitting every model
# model_names = nbd_mcmc.model_names

# Otherwise, enumerate the models to fit here
num_confs_list = [3]
"""The numbers of conformations to attempt to fit to the data."""

nsteps = 500000
"""The number of steps in each chain."""

num_temps = 8
"""The number of temperatures to run. Actual cores used will be num_temps + 1"""

num_chains = 10
"""The number of chains to run for each model."""

queue = 'parallel'
"""The LSF queue to submit the jobs to."""

time_limit = '24:00'
"""The estimated runtime of the job."""

#nbd_sites = ['c120', 'c122', 'c126', 'c15', 'c175', 'c179', 'c188', 'c36',
# 'c40', 'c47', 'c5', 'c54', 'c62', 'c68', 'c79']
nbd_sites = ['c5', 'c15', 'c47', 'c54', 'c62', 'c68', 'c79', 'c120', 'c122',
             'c126', 'c175', 'c179', 'c188']
"""The NBD sites to attempt to fit."""

num_replicates = 4
"""The number of replicates for each NBD mutant."""

python_str = '/home/jab69/virtualenvs/pysb/bin/python'

script_str = '/home/jab69/projects/tBidBaxLipo/tbidbaxlipo/mcmc/nbd_plate_mcmc_mpi_run.py'

if __name__ == '__main__':
    cmd_list = []
    for nbd_site in nbd_sites:
        for replicate in range(num_replicates):
            for num_confs in num_confs_list:
                for i in range(0, num_chains):
                    output_filename = '%s_rep%d_%dconfs_%d_s%d.out' % \
                                     (nbd_site, replicate, num_confs, nsteps, i)
                    cmd_list = ['bsub', '-a', 'openmpi',
                                '-n', str(num_temps+1),
                                '-W', time_limit,
                                '-q', queue,
                                '-o', output_filename,
                                'mpirun.lsf',
                                python_str,
                                script_str,
                                'random_seed=%d' % i,
                                "num_confs=%d" % num_confs,
                                "nbd_site=%s" % nbd_site,
                                "replicate=%d" % replicate,
                                "nsteps=%d" % nsteps]
                    print ' '.join(cmd_list)
                    subprocess.call(cmd_list)


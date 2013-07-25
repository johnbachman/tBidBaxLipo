"""
A script for submitting multiple MCMC jobs for the NBD plate data to LSF using
MPI-based parallel tempering.  Allows fitting models with alternative numbers
of assumed conformational states to different replicates for each NBD-labeled
Bax mutant.
"""

from tbidbaxlipo.mcmc import nbd_plate_mcmc
import subprocess
import numpy as np
import itertools

# == Models ==
num_confs_list = [2, 3, 4]
"""The numbers of conformations to fit for multiconf models."""

num_exps_list = [1, 2, 3]
"""The number of exponentials to fit for exponential models."""

multiconf_arg_list = ['model=multiconf num_confs=%d' % num_confs
                      for num_confs in num_confs_list]
exponential_arg_list = ['model=exponential num_exponentials=%d' % num_exps
                        for num_exps in num_exps_list]
model_arg_list = multiconf_arg_list + exponential_arg_list

# == Data ==
dataset = 'pti'
"""The NBD dataset to draw from."""

nbd_sites = ['c3']
"""The NBD sites to attempt to fit."""
#nbd_sites = ['c120', 'c122', 'c126', 'c15', 'c175', 'c179', 'c188', 'c36',
# 'c40', 'c47', 'c5', 'c54', 'c62', 'c68', 'c79']
#nbd_sites = ['c5', 'c15', 'c47', 'c54', 'c62', 'c68', 'c79', 'c120', 'c122',
#             'c126', 'c175', 'c179', 'c188']
nbd_site_arg_list = ['nbd_site=%s' % site for site in nbd_sites]

num_replicates = 3
"""The number of replicates for each NBD mutant."""
replicate_arg_list = ['replicate=%d' % i for i in range(num_replicates)]

# == Fitting ==
nsteps = 500000
"""The number of steps in each chain."""

num_temps = 8
"""The number of temperatures to run. Actual cores used will be num_temps + 1"""

num_chains = 10
"""The number of chains to run for each model."""
random_seed_arg_list = ['random_seed=%d' % i for i in range(num_chains)]

queue = 'parallel'
"""The LSF queue to submit the jobs to."""

time_limit = '24:00'
"""The estimated runtime of the job."""

python_str = '/home/jab69/virtualenvs/pysb/bin/python'
#python_str = 'python'
"""Path to the python executable."""

script_str = '/home/jab69/projects/tBidBaxLipo/tbidbaxlipo/mcmc/' \
             'nbd_plate_mcmc_mpi_run.py'
#script_str = 'script.py'
"""Path to the fitting script."""

def base_cmd_list(output_filename):
    """Get the base command list with the given output filename."""
    base_cmd_list = ['bsub', '-a', 'openmpi',
                '-n', str(num_temps+1),
                '-W', time_limit,
                '-q', queue,
                '-o', output_filename,
                'mpirun.lsf',
                python_str,
                script_str]
    return base_cmd_list

def output_filename_from_args(args):
    """Get the appropriate output filename given the current args."""
    # Join and then re-split the list at the spaces
    # This makes the string 'model=%s num_xxx=%d' into two separate args
    arg_strings = ' '.join(args).split(' ')
    # Now build up the list of key/val pairs and make a dict
    arg_dict = dict(arg_string.split('=') for arg_string in arg_strings)
    # Build and return the output filename
    if arg_dict['model'] == 'exponential':
        return '%s_%s_rep%s_%sexp_%d_s%s.out' % \
               (dataset, arg_dict['nbd_site'], arg_dict['replicate'],
                arg_dict['num_exponentials'], nsteps, arg_dict['random_seed'])
    elif arg_dict['model'] == 'multiconf':
        return '%s_%s_rep%s_%sconfs_%d_s%s.out' % \
               (dataset, arg_dict['nbd_site'], arg_dict['replicate'],
                arg_dict['num_confs'], nsteps, arg_dict['random_seed'])
    else:
        raise Exception("Invalid keyword argument for model.")

# Iterate over the Cartesian product of the different argument lists
for args in itertools.product(model_arg_list, random_seed_arg_list,
                              replicate_arg_list, nbd_site_arg_list):
    cmd_list = base_cmd_list(output_filename_from_args(args)) + list(args)
    print ' '.join(cmd_list)
    subprocess.call(cmd_list)


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
models = ['bax_heat',
          'bax_heat_reversible',
          'bax_heat_dimer',
          'bax_heat_dimer_reversible',
          'bax_heat_auto',
          'bax_heat_auto_reversible_activation',
          'bax_heat_auto_reversible',
          'bax_heat_auto_dimer',
          'bax_heat_auto_dimer_reversible',
          'bax_schwarz',
          'bax_schwarz_reversible',
          'bax_schwarz_dimer',
          'bax_schwarz_dimer_reversible',
          ]

cpt_types = ['one_cpt', 'lipo_sites']

model_arg_list = ['model=%s' % model_name for model_name in models]
cpt_type_arg_list = ['cpt_type=%s' % cpt_type for cpt_type in cpt_types]

nsteps = 50000
"""The number of steps in each chain."""

#num_temps = 8
#"""The number of temperatures to run. Actual cores used will be num_temps + 1"""

num_chains = 10
"""The number of chains to run for each model."""

random_seed_arg_list = ['random_seed=%d' % i for i in range(num_chains)]

queue = 'short'
"""The LSF queue to submit the jobs to."""

time_limit = '12:00'
"""The estimated runtime of the job."""

python_str = '/home/jab69/virtualenvs/pysb/bin/python'
"""Path to the python executable."""

script_str = '/home/jab69/projects/tBidBaxLipo/tbidbaxlipo/mcmc/' \
             'pore_mcmc.py'
"""Path to the fitting script."""

dataset_name = '130614'

def base_cmd_list(output_filename):
    """Get the base command list with the given output filename."""
    base_cmd_list = ['bsub', 
                '-W', time_limit,
                '-q', queue,
                '-o', output_filename,
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
    return '%s_%s_%s_%d_s%s.out' % \
           (dataset_name,
            arg_dict['cpt_type'],
            arg_dict['model'],
            nsteps,
            arg_dict['random_seed'])

# Iterate over the Cartesian product of the different argument lists
for args in itertools.product(model_arg_list,
                              cpt_type_arg_list,
                              random_seed_arg_list):
    fixed_args = ['nsteps=%d' % nsteps]
    cmd_list = base_cmd_list(output_filename_from_args(args)) + \
               list(args) + fixed_args
    print ' '.join(cmd_list)
    subprocess.call(cmd_list)

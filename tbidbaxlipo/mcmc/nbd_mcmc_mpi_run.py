"""
Example template for how to implement parallel tempering algorithm using MPI.

Run with, e.g.::

    mpiexec -np 5 python hello_mpi.py
"""

import sys
import numpy as np
from mpi4py import MPI

from bayessb import MCMCOpts
from bayessb.mpi.pt_mpi import PT_MPI_Master, PT_MPI_Worker
import tbidbaxlipo.plots.nbd_analysis as nbd
from tbidbaxlipo.mcmc.nbd_mcmc import model_names, nbd_site_names, NBD_MCMC
from pysb.integrate import Solver

if __name__ == '__main__':
    # Set the type of model to be built here
    from tbidbaxlipo.models.one_cpt import Builder

    # Prepare the data
    # ================

    tspan = nbd.time_other
    nbd_avgs, nbd_stds = nbd.calc_avg_std(normalize=True)

    # The c62 data is on a slightly different timescale and has two extra
    # datapoints. For simplicity, we simply remove them.
    nbd_avgs[1] = nbd_avgs[1][0:-2]
    nbd_stds[1] = nbd_stds[1][0:-2]

    # Parse the args
    # ==============
    # Keyword args are set at the command line as e.g., key=val
    # and subsequently split at the equals sign
    kwargs = dict([arg.split('=') for arg in sys.argv[1:]])

    # Arguments that can take a list of values (e.g., the nbd sites or
    # nbd observables) are expected to be separated by this character:
    SEP_CHAR = '-'

    # We set these all to None so later on we can make sure they were
    # properly initialized.
    nbd_sites = None
    nbd_observables = None
    random_seed = None
    model = None

    print "Keyword arguments: "
    print kwargs

    # Before we begin, we make sure we have all the keyword arguments that
    # we are going to need.
    if 'nbd_sites' not in kwargs or \
            'random_seed' not in kwargs or \
            'nsteps' not in kwargs or \
            'model' not in kwargs or \
            'nbd_observables' not in kwargs:
        raise Exception('One or more needed arguments was not specified! ' \
                'Arguments must include nbd_sites, random_seed, model, ' \
                'nbd_observables and nsteps.')

    # Because the NBD site(s) to fit has to be specified when the model
    # builder object is created, we get this arg first.
    # If more than one NBD site has been specified, it should be separated:
    nbd_sites = kwargs['nbd_sites'].split(SEP_CHAR)
    for site in nbd_sites:
        if site not in nbd_site_names:
            raise Exception('%s is not an allowed NBD site!' % site)

    # Now that we have the NBD site we can instantiate the Builder:
    builder = Builder(nbd_sites=nbd_sites)

    # The observable associated with the NBD site signal also has to be
    # specified:
    observables = [o.name for o in builder.model.observables]

    # Again, parse in case there is more than one observable:
    nbd_observables = kwargs['nbd_observables'].split(SEP_CHAR)
    for nbd_obs in nbd_observables:
        if nbd_obs not in observables:
            raise Exception('%s is not an allowed NBD observable!' % nbd_obs)

    # Make sure the number of nbd sites matches the number of observables:
    if len(nbd_sites) != len(nbd_observables):
        raise Exception('The number of nbd_sites must match the number of '
                        'nbd_observables!')

    # ...and then we get the model, which is specified as a string from the
    # set seen below.
    if kwargs['model'] not in model_names:
        raise Exception('%s is not an allowed model!' %
                        kwargs['model'])
    else:
        # Here we use a bit of Python trickery to avoid a long list of
        # if/elif statements: we append the model abbreviation to the
        # build_model function name and then eval it:
        build_model_cmd = 'builder.build_model_%s()' % kwargs['model']
        eval(build_model_cmd)
        model = builder.model

    # Set the random seed:
    random_seed = int(kwargs['random_seed'])

    # Set the number of steps:
    nsteps = int(kwargs['nsteps'])

    # A sanity check to make sure everything worked:
    if None in [model, nbd_sites, random_seed, nsteps, nbd_observables]:
        raise Exception('Something went wrong! One of the arguments to ' \
                        'do_fit was not initialized properly.')

    # We set the random_seed here because it affects our choice of initial
    # values
    np.random.seed(random_seed)

    # -----------------------------------------------------------------------
    # The communicator to use
    comm = MPI.COMM_WORLD
    # Number of chains/workers in the whole pool
    num_chains = comm.Get_size()
    # The rank of this chain (0 is the master, others are workers)
    rank = comm.Get_rank()
    # Forces the solver to use inline without testing first
    Solver._use_inline = True

    # Frequency for proposing swaps
    swap_period = 5
    # Temperature range
    min_temp = 1
    max_temp = 1e5

    # Create temperature array based on number of workers (excluding master)
    temps = np.logspace(np.log10(min_temp), np.log10(max_temp), num_chains-1)

    rate_step_sizes = np.logspace(np.log10(5e-4), np.log10(5e-2), num_chains-1)
    scaling_step_sizes = np.logspace(np.log10(1e-4), np.log10(1e-2), num_chains-1)

    # Initialize the MCMC arguments
    opts = MCMCOpts()
    opts.model = model
    opts.tspan = nbd.time_other
    opts.estimate_params = builder.estimate_params
    opts.initial_values = builder.random_initial_values()
    opts.nsteps = nsteps

    #opts.norm_step_size = np.array([0.01] + \
    #                    ([0.05] * (len(opts.estimate_params)-1)))
    opts.norm_step_size = np.array([scaling_step_sizes[rank-1]] + \
                        ([rate_step_sizes[rank-1]] * (len(opts.estimate_params)-1)))
    print "Step size: %g" % (opts.norm_step_size[1])

    opts.sigma_step = 0 # Don't adjust step size
    opts.sigma_max = 50
    opts.sigma_min = 0.01
    opts.accept_rate_target = 0.23
    opts.accept_window = 100
    opts.sigma_adj_interval = nsteps # Don't adjust step size

    opts.anneal_length = 0 # necessary so cooling does not occur
    opts.use_hessian = False
    opts.hessian_scale = 1
    opts.hessian_period = opts.nsteps / 20 #10
    opts.seed = random_seed
    opts.T_init = temps[rank - 1] # Use the temperature for this worker
    opts.thermo_temp = 1
    mcmc = NBD_MCMC(opts, nbd_avgs, nbd_stds, nbd_sites, nbd_observables,
                    builder)
    mcmc.initialize()

    # The master coordinates when swaps occur ---------
    if rank == 0:
        pt = PT_MPI_Master(comm, rank, opts, swap_period, num_chains)
        pt.run()
    # Everyone else runs MCMC steps and swaps when told -----------
    else:
        pt = PT_MPI_Worker(comm, rank, mcmc, swap_period)
        pt.run()




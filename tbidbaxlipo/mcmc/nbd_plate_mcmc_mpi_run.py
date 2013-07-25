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
from tbidbaxlipo.mcmc.nbd_plate_mcmc import NBDPlateMCMC
from tbidbaxlipo.data.nbd_plate_data import data, nbd_names
from pysb.integrate import Solver

if __name__ == '__main__':
    # Parse the args
    # ==============
    # Keyword args are set at the command line as e.g., key=val
    # and subsequently split at the equals sign
    kwargs = dict([arg.split('=') for arg in sys.argv[1:]])

    print "Keyword arguments: "
    print kwargs

    # Before we begin, we make sure we have all the keyword arguments that
    # we are going to need.
    if 'random_seed' not in kwargs or \
       'nsteps' not in kwargs or \
       'nbd_site' not in kwargs or \
       'replicate' not in kwargs or \
       'dataset' not in kwargs or \
       'model' not in kwargs:
        raise Exception('One or more needed arguments was not specified! ' \
                'Arguments must include random_seed, nsteps, model, ' \
                'nbd_site, dataset and replicate.')

    # Set the random seed:
    random_seed = int(kwargs['random_seed'])

    # Set the number of steps:
    nsteps = int(kwargs['nsteps'])

    # Prepare the data
    # ================
    # Figure out which dataset we're supposed to use
    dataset = kwargs['dataset']
    if dataset == 'plate':
        data = nbd_plate_data.data
        nbd_names = nbd_plate_data.nbd_names
        dataset_name = 'plate_%s_rep%d' % (nbd_site, replicate)
    elif dataset == 'pti':
        data = nbd_data.data
        nbd_names = nbd_data.nbd_names
        dataset_name = 'pti_%s_rep%d' % (nbd_site, replicate)
    else:
        raise Exception('Allowable values for dataset: plate, pti.')

    # Get the NBD mutant for the data we want to fit
    nbd_site = kwargs['nbd_site']
    if nbd_site not in nbd_names:
        raise Exception('%s not an allowable nbd_site.' % nbd_site)

    # Get the replicate to fit
    replicate = int(kwargs['replicate'])

    # Choose which data/replicate to fit
    tc = data[(nbd_site, replicate)]
    time = tc[:, 'TIME'].values
    values = tc[:, 'VALUE'].values

    # Prepare the model
    # =================
    model = kwargs['model']
    if model not in ['multiconf', 'exponential']:
        raise Exception("Model must be one of: multiconf, exponential.")
    # Multiconf models
    if model == 'multiconf':
        if 'num_confs' not in kwargs:
            raise Exception("Argument num_confs must be specified for model"
                            " of type multiconf.")
        num_confs = int(kwargs['num_confs'])
        b = multiconf.Builder()
        b.build_model_multiconf(num_confs, values[0])
    # Multi-exponential models
    elif model == 'exponential':
        if 'num_exponentials' not in kwargs:
            raise Exception("Argument num_exponentials must be specified for "
                            "model of type exponential.")
        num_exponentials = int(kwargs['num_exponentials'])
        b = exponential.Builder()
        b.build_model_exponential(num_exponentials, values[0])
    # This should never happen
    else:
        assert False

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

    rate_step_sizes = np.logspace(np.log10(5e-4), np.log10(5e-2),
                                  num_chains-1)
    scaling_step_sizes = np.logspace(np.log10(1e-4), np.log10(1e-2),
                                     num_chains-1)

    # Initialize the MCMC arguments
    opts = MCMCOpts()
    opts.model = b.model
    opts.tspan = time
    opts.estimate_params = b.estimate_params
    opts.initial_values = b.random_initial_values()
    opts.nsteps = nsteps

    #opts.norm_step_size = np.array([0.01] + \
    #                    ([0.05] * (len(opts.estimate_params)-1)))
    # FIXME not correct ordering
    #opts.norm_step_size = np.array([scaling_step_sizes[rank-1]] + \
    #                    ([rate_step_sizes[rank-1]] *
    #                     (len(opts.estimate_params)-1)))
    #print "Step size: %g" % (opts.norm_step_size[1])
    opts.norm_step_size = 0.05

    opts.sigma_step = 0 # Don't adjust step size
    #opts.sigma_max = 50
    #opts.sigma_min = 0.01
    #opts.accept_rate_target = 0.23
    #opts.accept_window = 100
    #opts.sigma_adj_interval = nsteps # Don't adjust step size

    opts.anneal_length = 0 # necessary so cooling does not occur
    opts.use_hessian = True
    opts.hessian_scale = 1
    opts.hessian_period = opts.nsteps / 20 #10
    opts.seed = random_seed
    opts.T_init = temps[rank - 1] # Use the temperature for this worker
    #opts.thermo_temp = 1

    mcmc = NBDPlateMCMC(opts, values, dataset_name, b, num_confs)

    mcmc.initialize()

    # The master coordinates when swaps occur ---------
    if rank == 0:
        pt = PT_MPI_Master(comm, rank, opts, swap_period, num_chains)
        pt.run()
    # Everyone else runs MCMC steps and swaps when told -----------
    else:
        pt = PT_MPI_Worker(comm, rank, mcmc, swap_period)
        pt.run()




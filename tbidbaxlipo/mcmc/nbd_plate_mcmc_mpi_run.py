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
from tbidbaxlipo.mcmc.nbd_plate_mcmc import NBDPlateMCMC, \
                                            parse_command_line_args
from tbidbaxlipo.models.nbd import multiconf, exponential
from tbidbaxlipo.data import nbd_data, nbd_plate_data
from pysb.integrate import Solver

if __name__ == '__main__':
    args = parse_command_line_args(sys.argv)

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
    b = args['builder']
    opts = MCMCOpts()
    opts.model = b.model
    opts.tspan = args['time']
    opts.estimate_params = b.estimate_params
    opts.initial_values = b.random_initial_values()
    opts.nsteps = args['nsteps']

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
    opts.seed = args['random_seed']
    opts.T_init = temps[rank - 1] # Use the temperature for this worker
    #opts.thermo_temp = 1

    mcmc = NBDPlateMCMC(opts, args['values'], args['dataset_name'], b)

    mcmc.initialize()

    # The master coordinates when swaps occur ---------
    if rank == 0:
        pt = PT_MPI_Master(comm, rank, opts, swap_period, num_chains)
        pt.run()
    # Everyone else runs MCMC steps and swaps when told -----------
    else:
        pt = PT_MPI_Worker(comm, rank, mcmc, swap_period)
        pt.run()




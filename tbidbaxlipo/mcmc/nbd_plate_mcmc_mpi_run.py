"""
Example template for how to implement parallel tempering algorithm using MPI.

Run with, e.g.::

    mpiexec -np 5 python hello_mpi.py
"""

import sys
import numpy as np

from bayessb import MCMCOpts
from bayessb.mpi.pt_mpi import PT_MPI_Master, PT_MPI_Worker
from mpi4py import MPI
from tbidbaxlipo.mcmc.nbd_plate_mcmc import NBDPlateMCMC, \
                                            parse_command_line_args,
                                            get_mcmc_opts
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
    opts = get_mcmc_opts(b, args, T_init=temps[rank - 1])

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




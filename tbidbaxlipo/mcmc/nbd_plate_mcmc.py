from tbidbaxlipo.models.nbd import multiconf, exponential
from tbidbaxlipo.data import nbd_data, nbd_plate_data
import bayessb
import numpy as np
import tbidbaxlipo.mcmc
from matplotlib import pyplot as plt
import pickle
import sys
import math
from bayessb.mpi.pt_mpi import PT_MPI_Master, PT_MPI_Worker
from mpi4py import MPI
import subprocess
import itertools
from collections import OrderedDict

###############################################
# MCMC class                                  #
###############################################

class NBDPlateMCMC(tbidbaxlipo.mcmc.MCMC):
    """ Document me"""
    def __init__(self, options, data, dataset_name, builder):
        tbidbaxlipo.mcmc.MCMC.__init__(self, options, builder)

        # Store data/model info
        self.data = data
        self.dataset_name = dataset_name

        # Set the MCMC functions
        self.options.likelihood_fn = self.likelihood
        self.options.prior_fn = self.builder.prior
        self.options.step_fn = self.step

    @staticmethod
    def likelihood(mcmc, position):
        mcmc.simulate(position=position)
        nbd = mcmc.solver.yexpr['NBD']
        err = np.sum((mcmc.data - nbd)**2 / (2 * ((0.02 * mcmc.data)**2)))
        return err

    def plot_data(self, axis):
        axis.plot(self.options.tspan, self.data)

    def get_observable_timecourses(self, position):
        timecourses = {}
        self.simulate(position=position)
        predicted_nbd = self.solver.yexpr['NBD']
        timecourses['Predicted NBD Signal'] = [self.options.tspan,
                                               predicted_nbd]
        return timecourses

    def get_residuals(self, position):
        """Return the residuals for the fit at the given position."""
        timecourses = self.get_observable_timecourses(position)
        (time, values) = timecourses['Predicted NBD Signal']
        return [time, self.data - values]

    def get_basename(self):
        return '%s_%s_%d_T%.2f_s%d' % (self.dataset_name,
                               self.builder.model.name,
                               self.options.nsteps,
                               self.options.T_init,
                               self.options.seed)

###############################################
# Utility functions                           #
###############################################

def parse_command_line_args(argv):
    # Keyword args are set at the command line as e.g., key=val
    # and subsequently split at the equals sign
    kwargs = dict([arg.split('=') for arg in argv])

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
        print ('One or more needed arguments was not specified! ' \
                'Arguments must include random_seed, nsteps, model, ' \
                'nbd_site, dataset and replicate.')
        sys.exit()

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
    elif dataset == 'pti':
        data = nbd_data.data
        nbd_names = nbd_data.nbd_names
    else:
        raise Exception('Allowable values for dataset: plate, pti.')

    # Get the name of the NBD mutant for the data we want to fit, and
    # check that the specified mutant is in this dataset
    nbd_site = kwargs['nbd_site']
    if nbd_site not in nbd_names:
        raise Exception('%s not an allowable nbd_site for dataset %s.' % \
                        (nbd_site, dataset))

    # Get the replicate to fit
    replicate = int(kwargs['replicate'])

    # Set the dataset name (for use in filenames, etc.)
    dataset_name = 'pti_%s_rep%d' % (nbd_site, replicate)

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

    return {'builder': b, 'random_seed': random_seed,
            'time': time, 'values': values, 'nsteps': nsteps,
            'dataset_name': dataset_name}

def get_mcmc_opts(builder, args, T_init=1):
    """Fills out the fields of the MCMCOpts object."""
    opts = bayessb.MCMCOpts()
    opts.model = builder.model
    opts.tspan = args['time']
    opts.estimate_params = builder.estimate_params
    opts.initial_values = builder.random_initial_values()
    opts.nsteps = args['nsteps']
    opts.norm_step_size = 0.1
    opts.sigma_step = 0
    #opts.sigma_max = 50
    #opts.sigma_min = 0.01
    #opts.accept_rate_target = 0.23
    #opts.accept_window = 100
    #opts.sigma_adj_interval = 200
    opts.anneal_length = 0
    opts.use_hessian = True # CHECK
    opts.hessian_scale = 1
    opts.hessian_period = opts.nsteps / 20 #10
    opts.seed = args['random_seed']
    opts.T_init = T_init
    return opts

###############################################
# Run scripts                                 #
###############################################

def run_single(argv):
    args = parse_command_line_args(argv)

    b = args['builder']
    np.random.seed(args['random_seed'])
    opts = get_mcmc_opts(b, args)

    from tbidbaxlipo.mcmc.nbd_plate_mcmc import NBDPlateMCMC
    mcmc = NBDPlateMCMC(opts, args['values'], args['dataset_name'], b)

    mcmc.do_fit()

    # Pickle it
    output_file = open('%s.mcmc' % mcmc.get_basename(), 'w')
    pickle.dump(mcmc, output_file)
    output_file.close()

    print "Done."

def run_parallel(mcmc, argv):
    """Run a parallel tempering job."""
    args = parse_command_line_args(argv)
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

###############################################
# Submit scripts                              #
###############################################

def output_filename_from_args(args):
    """Get the appropriate output filename given the current args."""
    # Join and then re-split the list at the spaces
    # This makes the string 'model=%s num_xxx=%d' into two separate args
    arg_strings = ' '.join(args).split(' ')
    # Now build up the list of key/val pairs and make a dict
    arg_dict = OrderedDict(arg_string.split('=') for arg_string in arg_strings)
    # Build and return the output filename
    output_filename = '_'.join(arg_dict.values()) + '.out'
    return output_filename

def submit_single():
    """Submit multiple MCMC jobs for the NBD plate data to LSF.

    Allows fitting models with alternative numbers of assumed conformational
    states to different replicates for each NBD-labeled Bax mutant.
    """
    # The numbers of conformations to attempt to fit to the data.
    num_confs_list = [2, 3, 4, 5]
    num_confs_args = ['num_confs=%d' % num_confs for num_confs in num_confs_list]
    #The NBD sites to attempt to fit."""
    #nbd_sites = ['c120', 'c122', 'c126', 'c15', 'c175', 'c179', 'c188', 'c36',
    # 'c40', 'c47', 'c5', 'c54', 'c62', 'c68', 'c79']
    nbd_sites = ['c175', 'c179', 'c5', 'c15']
    nbd_site_args = ['nbd_site=%s' % nbd_site for nbd_site in nbd_sites]
    # The number of replicates for each NBD mutant."""
    num_replicates = 4
    replicate_args = ['replicate=%d' % rep for rep in range(num_replicates)]
    # The number of chains to run for each model.
    num_chains = 10
    chain_index_args = ['random_seed=%d' % i for i in range(num_chains)]

    # The number of steps in each chain.
    nsteps = 50000
    # The LSF queue to submit the jobs to.
    queue = 'mini'
    # The estimated runtime of the job.
    time_limit = '00:10'

    def base_cmd_list(output_filename):
        base_cmd_list = ['bsub',
                '-W', time_limit,
                '-q', queue,
                '-o', output_filename,
                'python',
                '-m',
                'tbidbaxlipo.mcmc.nbd_plate_mcmc']
        return base_cmd_list

    fixed_args = ['nsteps=%d' % nsteps]

    for var_args in itertools.product(num_confs_args, nbd_site_args, replicate_args,
                                      chain_index_args):
        all_args = list(var_args) + fixed_args
        cmd_list = base_cmd_list(output_filename_from_args(all_args)) + all_args
        print ' '.join(cmd_list)
        subprocess.call(cmd_list)


def submit_parallel():
    pass

###############################################
# Main                                        #
###############################################

def main():
    usage =  '\nUsage:\n\n'
    usage += 'python nbd_plate_mcmc.py run_single [args]\n'
    usage += '  Run a single MCMC chain with the args in [args].\n'
    usage += 'python nbd_plate_mcmc.py run_parallel [args]\n'
    usage += '  Run a parallel tempering MCMC with the args in [args].\n'
    usage += 'python nbd_plate_mcmc.py submit_single\n'
    usage += '  Submit a set of single-chain jobs on Orchestra.\n'
    usage += 'python nbd_plate_mcmc.py submit_parallel\n'
    usage += '  Submit a set of parallel tempering jobs on Orchestra.\n'

    if len(sys.argv) <= 1:
        print usage
        sys.exit()

    if sys.argv[1] == 'run_single':
        run_single(sys.argv[2:])
    elif sys.argv[1] == 'run_parallel':
        run_parallel(sys.argv[2:])
    elif sys.argv[1] == 'submit_single':
        submit_single()
    elif sys.argv[1] == 'submit_parallel':
        submit_parallel()
    else:
        print usage
        sys.exit()

if __name__ == '__main__':
    main()

###############################################
# Tests                                       #
###############################################

def get_NBDPlateMCMC_instance():
    # Choose which data/replicate to fit
    tc = nbd_plate_data.data[('c68', 0)]
    time = tc[:, 'TIME'].values
    values = tc[:, 'VALUE'].values

    # Choose which model to build
    num_confs = 2
    b = multiconf.Builder()
    b.build_model_multiconf(num_confs, values[0])

    # Set initial estimates for scaling parameters
    scaling_parameters = [p for p in b.model.parameters
                          if p.name.endswith('_scaling')]
    for p in scaling_parameters:
        p.value = np.max(values)

    opts = bayessb.MCMCOpts()
    opts.model = b.model
    opts.tspan = time
    opts.estimate_params = [p for p in b.model.parameters
                            if not p.name.endswith('_0')]
    opts.initial_values = [p.value for p in opts.estimate_params]
    opts.nsteps = 10
    opts.T_init = 1
    opts.anneal_length = 0
    opts.use_hessian = False
    opts.sigma_step = 0
    opts.norm_step_size = 0.1
    opts.seed = 1

    mcmc = NBDPlateMCMC(opts, values, 'c68rep0', b)
    return mcmc

def test_NBDPlateMCMC_init():
    npm = get_NBDPlateMCMC_instance()
    assert True

def test_plot_data():
    npm = get_NBDPlateMCMC_instance()
    fig = plt.figure()
    axis = fig.gca()
    npm.plot_data(axis)
    #plt.show()
    assert True

def test_plot_fit():
    npm = get_NBDPlateMCMC_instance()
    npm.initialize()
    fig = npm.fit_plotting_function(position=npm.initial_position)
    fig.savefig('test_plot_fit.png')

def test_get_basename():
    npm = get_NBDPlateMCMC_instance()
    npm.get_basename()
    assert True


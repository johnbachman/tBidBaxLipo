from tbidbaxlipo.models.nbd_multiconf import Builder
from tbidbaxlipo.data.nbd_plate_data import data, nbd_names
import bayessb
import numpy as np
import tbidbaxlipo.mcmc
from matplotlib import pyplot as plt
import pickle
import sys
import math

class NBDPlateMCMC(tbidbaxlipo.mcmc.MCMC):
    """ Document me"""
    def __init__(self, options, data, dataset_name, builder, num_confs):
        tbidbaxlipo.mcmc.MCMC.__init__(self, options, builder)

        # Store data/model info
        self.data = data
        self.dataset_name = dataset_name
        self.num_confs = num_confs

        # Set the MCMC functions
        self.options.likelihood_fn = self.likelihood
        self.options.prior_fn = self.builder.prior
        self.options.step_fn = self.step

    # Pickling function for this class
    # Differs from superclass in that it assigns to self.likelihood
    # rather then self.get_likelihood_function()
    def __setstate__(self, state):
        bayessb.MCMC.__setstate__(self, state)
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

# TESTS #####
    def get_basename(self):
        return '%s_%dconfs_%d_T%.2f_s%d' % (self.dataset_name,
                               self.num_confs,
                               self.options.nsteps,
                               self.options.T_init,
                               self.options.seed)

def get_NBDPlateMCMC_instance():
    # Choose which data/replicate to fit
    tc = data[('c68', 0)]
    time = tc[:, 'TIME'].values
    values = tc[:, 'VALUE'].values

    # Choose which model to build
    num_confs = 2
    b = Builder()
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

    mcmc = NBDPlateMCMC(opts, values, 'c68rep0', b, num_confs)
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

# MAIN ######
if __name__ == '__main__':
    # Set the type of model to be built here
    from tbidbaxlipo.models.nbd_multiconf import Builder

    # Parse the args
    # ==============
    # Keyword args are set at the command line as e.g., key=val
    # and subsequently split at the equals sign
    kwargs = dict([arg.split('=') for arg in sys.argv[1:]])

    # We set these all to None so later on we can make sure they were
    # properly initialized.
    random_seed = None
    num_confs = None
    nsteps = None
    nbd_site = None
    replicate = None

    print "Keyword arguments: "
    print kwargs

    # Before we begin, we make sure we have all the keyword arguments that
    # we are going to need.
    if 'random_seed' not in kwargs or \
       'nsteps' not in kwargs or \
       'num_confs' not in kwargs or \
       'nbd_site' not in kwargs or \
       'replicate' not in kwargs:
        raise Exception('One or more needed arguments was not specified! ' \
                'Arguments must include random_seed, nsteps, num_confs, ' \
                'nbd_site, and replicate.')

    # Set the random seed:
    random_seed = int(kwargs['random_seed'])

    # Set the number of steps:
    nsteps = int(kwargs['nsteps'])

    # Get the NBD mutant for the data we want to fit
    nbd_site = kwargs['nbd_site']
    if nbd_site not in nbd_names:
        raise Exception('%s not an allowable nbd_site.' % nbd_site)

    # Get the replicate to fit
    replicate = int(kwargs['replicate'])

    # Get the number of conformations
    num_confs = int(kwargs['num_confs'])

    # A sanity check to make sure everything worked:
    if None in [random_seed, nsteps, num_confs, nbd_site, replicate]:
        raise Exception('Something went wrong! One of the arguments to ' \
                        'do_fit was not initialized properly.')

    # Prepare the data and model
    # ==========================
    # Choose which data/replicate to fit
    tc = data[(nbd_site, replicate)]
    time = tc[:, 'TIME'].values
    values = tc[:, 'VALUE'].values

    # Build the model
    b = Builder()
    b.build_model_multiconf(num_confs, values[0])

    # Build the MCMCOpts
    # ==================
    # We set the random_seed here because it affects our choice of initial
    # values
    np.random.seed(random_seed)

    opts = bayessb.MCMCOpts()
    opts.model = b.model
    opts.tspan = time
    opts.estimate_params = b.estimate_params
    #opts.estimate_params = [p for p in b.model.parameters
    #                        if not p.name.endswith('_0')]
    #opts.initial_values = [p.value for p in b.estimate_params]
    opts.initial_values = b.random_initial_values()
    opts.nsteps = nsteps
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
    opts.seed = random_seed
    opts.T_init = 1
    dataset_name = '%s_rep%d' % (nbd_site, replicate)

    from tbidbaxlipo.mcmc.nbd_plate_mcmc import NBDPlateMCMC
    mcmc = NBDPlateMCMC(opts, values, dataset_name, b, num_confs)

    mcmc.do_fit()

    # Pickle it
    output_file = open('%s.mcmc' % mcmc.get_basename(), 'w')
    pickle.dump(mcmc, output_file)
    output_file.close()

    print "Done."



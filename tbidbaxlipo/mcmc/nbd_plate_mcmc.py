from tbidbaxlipo.models.nbd import multiconf, exponential
from tbidbaxlipo.data import nbd_data, nbd_plate_data
import bayessb
import numpy as np
import tbidbaxlipo.mcmc
from matplotlib import pyplot as plt
import pickle
import sys
import math

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

# TESTS #####
    def get_basename(self):
        return '%s_%s_%d_T%.2f_s%d' % (self.dataset_name,
                               self.builder.model.name,
                               self.options.nsteps,
                               self.options.T_init,
                               self.options.seed)

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

def parse_command_line_args(argv):
    # Keyword args are set at the command line as e.g., key=val
    # and subsequently split at the equals sign
    kwargs = dict([arg.split('=') for arg in argv[1:]])

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

# MAIN ######
if __name__ == '__main__':
    args = parse_command_line_args(sys.argv)

    # Build the MCMCOpts
    # ==================
    b = args['builder']

    # We set the random_seed here because it affects our choice of initial
    # values
    np.random.seed(args['random_seed'])

    opts = bayessb.MCMCOpts()
    opts.model = b.model
    opts.tspan = args['time']
    opts.estimate_params = b.estimate_params
    #opts.estimate_params = [p for p in b.model.parameters
    #                        if not p.name.endswith('_0')]
    #opts.initial_values = [p.value for p in b.estimate_params]
    opts.initial_values = b.random_initial_values()
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
    opts.T_init = 1

    from tbidbaxlipo.mcmc.nbd_plate_mcmc import NBDPlateMCMC
    mcmc = NBDPlateMCMC(opts, args['values'], args['dataset_name'], b)

    mcmc.do_fit()

    # Pickle it
    output_file = open('%s.mcmc' % mcmc.get_basename(), 'w')
    pickle.dump(mcmc, output_file)
    output_file.close()

    print "Done."


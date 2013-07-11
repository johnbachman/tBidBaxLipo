from tbidbaxlipo.models.nbd_multiconf import Builder
from tbidbaxlipo.data.nbd_plate_data import data
import bayessb
import numpy as np
import tbidbaxlipo.mcmc
from matplotlib import pyplot as plt

class NBDPlateMCMC(tbidbaxlipo.mcmc.MCMC):
    """ Document me"""
    def __init__(self, options, data, dataset_name, builder, num_confs):
        tbidbaxlipo.mcmc.MCMC.__init__(self, options, builder)

        # Store data/model info
        self.data = data
        self.dataset_name = dataset_name
        self.num_confs = num_confs

        # Set the MCMC functions
        self.options.likelihood_fn = self.get_likelihood_function()
        self.options.prior_fn = self.builder.prior
        self.options.step_fn = self.step

    def get_likelihood_function(self):
        # This is a shitty hack that relies on the known ordering of the
        # parameters from the way the model is created by the builder.
        # According to this order,
        # the rate parameters will have indices 1, 3, ... 
        # in the cur_params position array; the scaling parameters will have
        # indices 2, 4, ....
        observables = self.options.model.observables

        def likelihood(mcmc, position):
            x = mcmc.simulate(position=position, observables=True)
            aggregate_observable = 0
            cur_position = mcmc.cur_params(position=position)
            for i, obs in enumerate(observables):
                if i == 0:
                    continue
                scaling_parameter = cur_position[i*2]
                aggregate_observable += x[obs.name] * scaling_parameter
            err = ((self.data - aggregate_observable)**2 /
                   (2 * ((0.1 * self.data)**2)))
            return err
        return likelihood

    def plot_data(self, axis):
        axis.plot(self.options.tspan, self.data)

    def get_observable_timecourses(self, position):
        timecourses = {}
        x = self.simulate(position=position, observables=True)
        observables = self.options.model.observables
        aggregate_observable = 0
        cur_position = self.cur_params(position=position)
        for i, obs in enumerate(observables):
            if i == 0:
                continue
            scaling_parameter = cur_position[i*2]
            aggregate_observable += x[obs.name] * scaling_parameter
        timecourses['Predicted NBD Signal'] = [self.options.tspan,
                                               aggregate_observable]
        return timecourses

    def get_basename(self):
        return '%s_%dconfs_%d_s%d' % (self.dataset_name,
                               self.num_confs,
                               self.options.nsteps,
                               self.options.seed)

def get_NBDPlateMCMC_instance():
    # Choose which data/replicate to fit
    tc = data[('c68', 0)]
    time = tc[:, 'TIME'].values
    values = tc[:, 'VALUE'].values
    # Normalize values by subtracting initial value
    values = values - values[0]

    # Choose which model to build
    num_confs = 2
    b = Builder()
    b.build_model_multiconf(num_confs)

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

"""
A class for fitting the mechanistic pore formation models to dye release data.
"""

import tbidbaxlipo.mcmc
from tbidbaxlipo.models.one_cpt import Builder
import bayessb
import pickle
import numpy as np

class PoreMCMC(tbidbaxlipo.mcmc.MCMC):
    """Fit mechanistic tBid/Bax models to dye release titration data.

    Initialize parent tbidbaxlipo.mcmc.MCMC and then set additional fields.

    Parameters
    ----------
    options : MCMCOpts
        Options for MCMC initialization.
    data :
        pandas data structure containing the timecourses for each Bax condition.
    builder : tbidbaxlipo.models.core.Builder
    """
    def __init__(self, options, data, builder):
        # Call the superclass constructor
        tbidbaxlipo.mcmc.MCMC.__init__(self, options, builder)

        # Set the MCMC functions
        self.options.likelihood_fn = self.get_likelihood_function()
        self.options.prior_fn = self.builder.prior
        self.options.step_fn = self.step

    # Implementations of necessary functions
    def get_likelihood_function(self):
        return None

    def plot_data(self, axis):
        raise NotImplementedError()

    def get_observable_timecourses(self, position):
        raise NotImplementedError()

    def get_basename(self):
        return '%s_%s_%d_s%d' % (self.dataset_name,
                                 self.options.model.name,
                                 self.options.nsteps,
                                 self.options.seed)

def get_PoreMCMC_instance():
    b = Builder()
    b.build_model_tar()

    data = pickle.load(open('../data/130614_norm_timecourses_df.pck'))

    opts = bayessb.MCMCOpts()
    opts.model = b.model
    opts.tspan = np.linspace(0, 6000, 196)
    opts.estimate_params = [p for p in b.model.parameters]
    opts.initial_values = [p.value for p in opts.estimate_params]
    opts.nsteps = 10
    opts.T_init = 1
    opts.anneal_length = 0
    opts.use_hessian = False
    opts.sigma_step = 0
    opts.norm_step_size = 0.1
    opts.seed = 1

    p = PoreMCMC(opts, data, b)

def test_PoreMCMC_init():
    get_PoreMCMC_instance()
    assert True

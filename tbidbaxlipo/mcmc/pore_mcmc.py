"""
A class for fitting the mechanistic pore formation models to dye release data.
"""

import tbidbaxlipo.mcmc
from tbidbaxlipo.models.one_cpt import Builder
import bayessb
import pickle
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import collections

model_names = ['bax_heat',
               'bax_heat_reversible',
               'bax_heat_dimer',
               'bax_heat_dimer_reversible',
               'bax_heat_auto',
               'bax_heat_auto_reversible',
               'bax_heat_auto_dimer',
               'bax_heat_auto_dimer_reversible',
              ]

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
    def __init__(self, options, data, dataset_name, builder):
        # Call the superclass constructor
        tbidbaxlipo.mcmc.MCMC.__init__(self, options, builder)

        # Store the data
        self.data = data
        self.dataset_name = dataset_name

        # Set the MCMC functions
        self.options.likelihood_fn = self.get_likelihood_function()
        self.options.prior_fn = self.builder.prior
        self.options.step_fn = self.step

    # Implementations of necessary functions
    def get_likelihood_function(self):

        def likelihood(mcmc, position):
            err = 0
            for bax_conc in mcmc.data.columns:
                # Get the data for this concentration
                tc = mcmc.data[bax_conc]
                y_data  = np.array(tc[:,'MEAN'])
                time = np.array(tc[:,'TIME'])
                mcmc.solver.tspan = time # set the time span

                # Get the simulated data for this concentration
                mcmc.options.model.parameters['Bax_0'].value = bax_conc
                x = mcmc.simulate(position=position, observables=True)
                avg_pores = x['pores']/ \
                            mcmc.options.model.parameters['Vesicles_0'].value
                y_mod = 1 - np.exp(-avg_pores)

                # Calculate the error, accounting for the SD at this
                # concentration.
                # Skip the first timepoint--the SD is 0 (due to normalization)
                # and hence gives nan when calculating the error.
                err += np.sum(((y_data[1:] - y_mod[1:])**2) / \
                        (2 * (np.array(tc[:,'SD'][1:]) ** 2)))
            return err
        return likelihood

    def plot_data(self, axis):
        # Plot the titration of Bax timecourses
        for bax_conc in self.data.columns:
            tc = self.data[bax_conc]
            axis.plot(tc[:,'TIME'], tc[:,'MEAN'], # error=tc[:,'SD'],
                       color='gray')

    def get_observable_timecourses(self, position):
        """Return the timecourses for all concentrations."""
        timecourses = collections.OrderedDict()

        for bax_conc in self.data.columns:
            # Get the timepoints for this concentration
            tc = self.data[bax_conc]
            time = np.array(tc[:,'TIME'])
            self.solver.tspan = time # set the time span
            self.options.model.parameters['Bax_0'].value = bax_conc
            x = self.simulate(position=position, observables=True)
            avg_pores = x['pores'] / \
                        self.options.model.parameters['Vesicles_0'].value
            y_mod = 1 - np.exp(-avg_pores)
            timecourses['Bax %d nM' % bax_conc] = [time, y_mod]
        return timecourses

    def get_basename(self):
        return '%s_%s_%s_%d_s%d' % (self.dataset_name,
                                 self.builder.get_module(),
                                 self.options.model.name,
                                 self.options.nsteps,
                                 self.options.seed)

def get_PoreMCMC_instance():
    b = Builder()
    b.build_model_bax_heat()

    data = pickle.load(open('../data/130614_norm_timecourses_df.pck'))

    opts = bayessb.MCMCOpts()
    opts.model = b.model
    opts.tspan = np.linspace(0, 6000, 196)
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

    pm = PoreMCMC(opts, data, dataset_name, b)
    return pm

def test_PoreMCMC_init():
    pm = get_PoreMCMC_instance()
    assert True

def test_plot_data():
    pm = get_PoreMCMC_instance()
    fig = plt.figure()
    axis = fig.gca()
    pm.plot_data(axis)
    #plt.show()
    assert True

def test_plot_fit():
    pm = get_PoreMCMC_instance()
    pm.initialize()
    fig = pm.fit_plotting_function(position=pm.initial_position)
    fig.savefig('test_plot_fit.png')
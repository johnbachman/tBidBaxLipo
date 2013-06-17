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

        # Store the data
        self.data = data

        # Set the MCMC functions
        self.options.likelihood_fn = self.likelihood
        self.options.prior_fn = self.builder.prior
        self.options.step_fn = self.step

    # Implementations of necessary functions
    def likelihood(self, position):
        for bax_conc in self.data_columns:
            # Get the data for this concentration
            tc = self.data[bax_conc]
            y_data  = np.array(tc[:,'MEAN'])
            time = np.array(tc[:,'TIME'])
            self.solver.tspan = time # set the time span

            # Get the simulated data for this concentration
            self.options.model.parameters['Bax_0'].value = bax_conc
            x = self.simulate(position=position, observables=True)
            avg_pores = x['pores']/ \
                        self.options.model.parameters['Vesicles_0'].value
            y_mod = 1 - np.exp(-avg_pores)

            # Calculate the error, accounting for the SD at this concentration
            err += np.sum((y_data - y_mod)**2) / \
                          (2 * (np.array(tc[:,'SD']) ** 2))

        return err

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

    def fit_plotting_function(self, position):
        """Gets the observable timecourse and plots it against the data."""
        # Run the simulation at the given position
        timecourses = self.get_observable_timecourses(position)
        # Make the plot
        fig = Figure()
        ax = fig.gca()
        # Add the data to the plot
        self.plot_data(ax)

        # Add the simulations to the plot
        for obs_name, timecourse in timecourses.iteritems():
            ax.plot(timecourse[0], timecourse[1], label=obs_name)

        # Label the plot
        ax.set_xlabel('Time')
        ax.set_ylabel('Concentration')
        fontP = FontProperties() 
        fontP.set_size('small')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='upper left', prop=fontP, ncol=1, bbox_to_anchor=(1, 1),
             fancybox=True, shadow=True)
        #ax.legend(loc='upper center', prop=fontP, ncol=5,
        #            bbox_to_anchor=(0.5, 1.1), fancybox=True, shadow=True)
        canvas = FigureCanvasAgg(fig)
        fig.set_canvas(canvas)
        return fig


    def get_basename(self):
        return '%s_%s_%d_s%d' % (self.dataset_name,
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

    pm = PoreMCMC(opts, data, b)
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

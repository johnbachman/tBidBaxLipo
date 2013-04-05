"""
Functions for fitting the mechanistic NBD insertion models to the NBD-Bax mutant
fluorescence data using MCMC.

.. todo:: Likelihood should be scaled by number of timepoints, otherwise it

overwhelms the prior for "no reason".

.. todo:: Ideally, would have a way of pre-equilibrating the system for the

just-Bax condition, and then perturb it with the addition of tBid.
"""

import bayessb
import os
from pysb.integrate import odesolve
import numpy
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
from tbidbaxlipo.util.report import Report
from matplotlib.font_manager import FontProperties
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

model_names = ['ta', 'tai', 'taid', 'taidt', 'tair', 'taird', 'tairdt',
               'tad', 'tadt', 'tar', 'tard', 'tardt']

nbd_site_names = ['c3', 'c62']

class NBD_MCMC(bayessb.MCMC):
    """Fit mechanistic tBid/Bax models to NBD data.

    Initialize parent bayessb.MCMC and then set additional fields.

    Parameters
    ----------
    options : MCMCOpts
        Options for MCMC initialization.
    nbd_avgs : numpy.array
        data mean
    nbd_stds : numpy.array
        data SD
    nbd_sites : list of strings
        Sites from data to fit.
    nbd_observables : list of strings
        Observables from model to fit to the sites in nbd_sites.
    builder : tbidbaxlipo.models.core.Builder
    """

    def __init__(self, options, nbd_avgs, nbd_stds, nbd_sites, nbd_observables,
                 builder):
        bayessb.MCMC.__init__(self, options)
        self.nbd_avgs = nbd_avgs
        self.nbd_stds = nbd_stds
        self.nbd_sites = nbd_sites
        self.nbd_observables = nbd_observables
        self.builder = builder

        # Set the MCMC functions
        self.options.likelihood_fn = self.get_likelihood_function()
        self.options.prior_fn = self.builder.prior
        self.options.step_fn = self.step

    # Pickling functions for this class
    def __getstate__(self):
        self.options.likelihood_fn = None
        self.options.prior_fn = None
        self.options.step_fn = None
        mcmc_state = bayessb.MCMC.__getstate__(self)
        #nbd_mcmc_state = self.__dict__.copy()
        #return (mcmc_state, nbd_mcmc_state)
        return mcmc_state

    def __setstate__(self, state):
        #(mcmc_state, nbd_mcmc_state) = state
        bayessb.MCMC.__setstate__(self, state)
        #self.__dict__.update(nbd_mcmc_state)
        self.options.likelihood_fn = self.get_likelihood_function()
        self.options.prior_fn = self.builder.prior
        self.options.step_fn = self.step

    # MCMC Functions
    # ==============
    def do_fit(self):
        """Runs MCMC on the given model."""

        self.initialize()

        # Print initial parameter values
        init_vals = zip([p.name for p in self.options.model.parameters],
                          self.cur_params(position=self.initial_position))
        init_vals_str = 'Initial values:\n'
        init_vals_str += '\n'.join(['%s: %g' % (init_vals[i][0],
                                                 init_vals[i][1])
                                     for i in range(0, len(init_vals))])
        print "------------------------"
        print init_vals_str
        print "------------------------"

        # Run it!
        self.run()

        # Pickle it, setting functions to None:
        #self.options.likelihood_fn = None
        #self.options.prior_fn = None
        #self.options.step_fn = None
        # Restore the functions for interactive use
        #self.options.likelihood_fn = likelihood
        #self.options.prior_fn = prior
        #self.options.step_fn = step

    def generate_figures(self, report_name='report', do_report=True,
                         mixed_start=None, num_samples=500):
        """Plots a series of useful visualizations of the walk, the
        quality of the fit, etc."""
        plt.ion()

        if do_report:
            rep = Report()

        # Plot "Before" curves -------
        before_fig = self.fit_plotting_function(self.initial_position)
        if do_report:
            rep.add_figure(before_fig)

        # Print initial parameter values
        init_vals = zip([p.name for p in self.options.model.parameters],
                          self.cur_params(position=self.initial_position))
        init_vals_str = 'Initial values:\n'
        init_vals_str += '\n'.join(['%s: %g' % (init_vals[i][0],
                                                init_vals[i][1])
                                     for i in range(0, len(init_vals))])
        print init_vals_str
        if do_report:
            rep.add_text(init_vals_str)

        # Plot "After" curves ------------
        # Set to last fit position
        if mixed_start is None:
            mixed_start = self.options.nsteps / 2
        mixed_positions = self.positions[mixed_start:,:]
        mixed_accepted_positions = mixed_positions[self.accepts[mixed_start:]]
        last_position = mixed_accepted_positions[-1,:]

        after_fig = self.fit_plotting_function(last_position)
        if do_report:
            rep.add_figure(after_fig)

        # Print final parameter values
        last_fit_params = self.cur_params(position=last_position)
        last_vals = zip([p.name for p in self.options.model.parameters],
                           last_fit_params)
        last_vals_str = 'Final values:\n'
        last_vals_str += '\n'.join(['%s: %g' % (last_vals[i][0],
                                                last_vals[i][1])
                                    for i in range(0, len(last_vals))])
        print last_vals_str
        if do_report:
            rep.add_text(last_vals_str)

        """
        # Plot sampling of fits # FIXME should be real sampling------
        plt.figure()
        plot_data(nbd_site)
        max_position_index = len(mixed_accepted_positions) - 1
        num_to_plot = min(num_samples, max_position_index)
        for i in range(num_to_plot):
            cur_position = mixed_accepted_positions[max_position_index - i,:]
            x = nbd_timecourse(mcmc, cur_position, nbd_observable)
            plt.plot(mcmc.options.tspan, x, 'b', label=nbd_site, alpha=0.02)
        plt.title("Last %d accepted positions" % num_to_plot)
        plt.show()
        if do_report:
            rep.add_current_figure()

        # Plot marginal distributions of parameters ---------------
        for i, cur_param in enumerate(mcmc.options.estimate_params):
            plt.figure()
            plt.hist(mixed_accepted_positions[:,i])
            plt.title("%s, last %d accepts: initval %f" %
                    (cur_param.name, len(mixed_accepted_positions[:,i]),
                     cur_param.value))
            plt.show()

            if do_report:
                rep.add_current_figure()

        # Plot convergence traces of all parameters
        plt.figure()
        plt.plot(mcmc.positions)
        plt.title("Parameter traces")
        plt.legend([p.name for p in mcmc.options.estimate_params],
                    loc='lower left', prop={'size':7})
        plt.show()
        if do_report:
            rep.add_current_figure()
        """

        # Add code to report
        if do_report:
            # Add the code for the fitting (this file)
            rep.add_python_code('nbd_mcmc_pysb.py')

            # Add the code for the model
            rep.add_python_code('models/core.py')
            rep.add_python_code('models/one_cpt.py')

            # Write the report 
            rep.write_report(report_name)

    def plot_data(self, axis):
        """Plots the current data into the given axis."""
        alpha = 0.5
        if 'c3' in self.nbd_sites:
            axis.plot(self.options.tspan, self.nbd_avgs[0], 'r.',
                     label='c3 data', alpha=alpha)
        if 'c62' in self.nbd_sites:
            axis.plot(self.options.tspan, self.nbd_avgs[1], 'g.',
                     label='c62 data', alpha=alpha)
        #plt.plot(nbd.time_other, nbd_avgs[2], 'b.', label='c120 data',
        #          alpha=alpha)
        #plt.plot(nbd.time_other, nbd_avgs[3], 'm.', label='c122 data',
        #         alpha=alpha)
        #plt.plot(nbd.time_other, nbd_avgs[4], 'k.', label='c126 data',
        #         alpha=alpha)

    def get_observable_timecourses(self, position):
        """Gets the timecourses associated with the experimental observables."""
        # TODO Need to be able to get the indices for scaling parameters
        # from the model so that they're not hardcoded

        # Run the simulation and scale the simulated timecourses
        yout = self.simulate(position, observables=True)
        params = self.cur_params(position)
        timecourses = []
        for obs in self.nbd_observables:
            timecourses.append((yout[obs] /
                           self.options.model.parameters['Bax_0'].value)
                          * params[3])

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
        for timecourse in timecourses:
            ax.plot(self.options.tspan, timecourse)
        # Label the plot
        ax.set_xlabel('Time')
        ax.set_ylabel('Concentration')
        fontP = FontProperties() 
        fontP.set_size('small')
        ax.legend(loc='upper center', prop=fontP, ncol=5,
                    bbox_to_anchor=(0.5, 1.1), fancybox=True, shadow=True)
        canvas = FigureCanvasAgg(fig)
        fig.set_canvas(canvas)
        return fig

    # A function to generate the likelihood function
    def get_likelihood_function(self):
        """Returns a likelihood function for the specified NBD site."""
        data_indices = []
        if 'c3' in self.nbd_sites:
            data_indices.append(0)
        if 'c62' in self.nbd_sites:
            data_indices.append(1)
        # Make sure the list is not empty
        if not data_indices:
            raise Exception('Failed to initialize data_indices!')

        # Make sure that there is a corresponding data trajectory for each
        # observable
        if len(data_indices) != len(self.nbd_observables):
            raise Exception('Length of list of nbd_sites does not match the '
                            'list of nbd_observables!')

        def likelihood(mcmc, position):
            yout = mcmc.simulate(position, observables=True)
            # TODO Need to be able to get the indices from the model so that
            # they're not hardcoded
            params = mcmc.cur_params(position)
            err = 0
            for data_index, obs_name in zip(data_indices, self.nbd_observables):
                timecourse = ((yout[obs_name] /
                               mcmc.options.model.parameters['Bax_0'].value)
                              * params[3])
                err += numpy.sum((self.nbd_avgs[data_index] - timecourse)**2 /
                                 (2 * self.nbd_stds[data_index]**2))
            return err

        return likelihood

    def get_basename(self):
        """A function for standardizing, in one place, the format for pickled
        NBD_MCMC objects.
        """
        return '%s_%s_%s_%d_thermo%.2f_s%d' % (self.options.model.name,
                                    '-'.join(self.nbd_sites),
                                    '-'.join(self.nbd_observables),
                                    self.options.nsteps,
                                    #self.options.T_init,
                                    np.log10(self.options.thermo_temp),
                                    self.options.seed)

    @staticmethod
    def step(mcmc):
        """The function to call at every iteration. Currently just prints
        out a few progress indicators.
        """
        window = mcmc.options.accept_window

        local_acc = numpy.sum(mcmc.accepts[(mcmc.iter - window):mcmc.iter]) / \
                              float(window)

        if mcmc.iter % 20 == 0:
            print 'iter=%-5d  sigma=%-.3f  T=%-.3f  loc_acc=%-.3f  ' \
                  'glob_acc=%-.3f  lkl=%g  prior=%g  post=%g' % \
                  (mcmc.iter, mcmc.sig_value, mcmc.T,
                   local_acc,
                   mcmc.acceptance/(mcmc.iter+1.), mcmc.accept_likelihood,
                   mcmc.accept_prior, mcmc.accept_posterior)

# Chain handling helper function
# ==============================
def import_mcmc_groups(filenames):
    """Loads the chains into groups representing multiple runs of the same
    model.

    Assumes that the pickle filenames are structured as

        basename = '%s_%s_%s_%d_s%d.xxx' % (model, nbd_sites, nbd_observables,
                                        nsteps, random_seed)

    With the suffix ``.xxx`` separated by a dot and the seed coming last in
    the underscore-separated arguments.

    Parameters
    ----------
    filenames : list of strings
        List of strings representing the chain filenames to be sorted into
        groups, e.g., of the type returned by ``glob.glob()``.

    Returns
    -------
    dict of lists of MCMC filenames
        The keys in the dict are the filename prefixes that represent the
        arguments to the MCMC procedure (e.g., ``tard_c3_iBax_4000``. Each dict
        entry contains a list of MCMC filenames associated with those run
        conditions.
    """

    mcmc_groups = {}

    for filename in filenames:
        # Split off the suffix from the filename
        (prefix, suffix) = filename.split('.')
        # Separate the filename into the final argument identifying the
        # random seed, and everything that comes before it:
        (mcmc_args, seed) = prefix.rsplit('_', 1)
        mcmc_args = os.path.basename(mcmc_args)

        # Check to see if we've imported another chain of this type already.
        # If so, add the current chain to the list of chains for this group:
        if mcmc_args in mcmc_groups:
            mcmc_groups[mcmc_args].append(filename)
        # If not, create a new entry in the dict containing this MCMC
        else:
            mcmc_groups[mcmc_args] = [filename]

    return mcmc_groups

"""Functions for fitting the linear NBD insertion models to the NBD-Bax mutant
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

model_names = ['ta', 'tai', 'taid', 'taidt', 'tair', 'taird', 'tairdt',
               'tad', 'tadt', 'tar', 'tard', 'tardt']

nbd_site_names = ['c3', 'c62']

class NBD_MCMC(bayessb.MCMC):

    def __init__(self, options, nbd_avgs, nbd_stds, nbd_site, nbd_observable,
                 builder):
        """Initialize parent bayessb.MCMC and then set additional fields.

        Parameters
        ----------
        data_avgs : list of numpy.array
        data_stds : list of numpy.array
        nbd_site : string
        nbd_observable : string
        builder : tbidbaxlipo.models.core.Builder
        """

        bayessb.MCMC.__init__(self, options)
        self.nbd_avgs = nbd_avgs
        self.nbd_stds = nbd_stds
        self.nbd_site = nbd_site
        self.nbd_observable = nbd_observable
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

    def nbd_timecourse(mcmc, position, nbd_observable):
        """Simulates the model at the given parameter position and returns
        the appropriately scaled timecourse for the given NBD site."""

        x = mcmc.simulate(position=position, observables=True)

        total_Bax = mcmc.options.model.parameters['Bax_0'].value
        cur_params = mcmc.cur_params(position=position)
        scaling_factor = cur_params[3]
        return (x[nbd_observable] / total_Bax) * scaling_factor 

    def generate_figures(mcmc, nbd_site, nbd_observable, do_report=True,
                         mixed_start=None, basename='nbd_mcmc',
                         num_samples=500):
        """Takes an MCMC chain and plots a series of useful visualizations of
        the walk, the quality of the fit, etc."""
        plt.ion()

        if do_report:
            rep = Report()

        # Plot "Before" curves -------
        x = nbd_timecourse(mcmc, mcmc.initial_position, nbd_observable)
        plt.figure()
        plot_data(nbd_site)
        plt.plot(mcmc.options.tspan, x, 'b', label=nbd_site)
        plt.legend(loc='lower right')
        plt.title('Before')
        plt.show()
        if do_report:
            rep.add_current_figure()

        # Print initial parameter values
        init_vals = zip([p.name for p in mcmc.options.model.parameters],
                          mcmc.cur_params(position=mcmc.initial_position))
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
            mixed_start = mcmc.options.nsteps / 2
        mixed_positions = mcmc.positions[mixed_start:,:]
        mixed_accepted_positions = mixed_positions[mcmc.accepts[mixed_start:]]
        last_position = mixed_accepted_positions[-1,:]

        x = nbd_timecourse(mcmc, last_position, nbd_observable)
        plt.figure()
        plot_data(nbd_site)
        plt.plot(mcmc.options.tspan, x, 'b', label=nbd_site)
        plt.legend(loc='lower right')
        plt.title('Final accepted position')
        plt.show()
        if do_report:
            rep.add_current_figure()

        # Print final parameter values
        last_fit_params = mcmc.cur_params(position=last_position)
        last_vals = zip([p.name for p in mcmc.options.model.parameters],
                           last_fit_params)
        last_vals_str = 'Final values:\n'
        last_vals_str += '\n'.join(['%s: %g' % (last_vals[i][0],
                                                last_vals[i][1])
                                    for i in range(0, len(last_vals))])
        print last_vals_str
        if do_report:
            rep.add_text(last_vals_str)

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

        # Add code to report
        if do_report:
            # Add the code for the fitting (this file)
            rep.add_python_code('nbd_mcmc_pysb.py')

            # Add the code for the model
            rep.add_python_code('models/core.py')
            rep.add_python_code('models/one_cpt.py')

            # Write the report 
            rep.write_report(basename)

    def plot_data(self):
        alpha = 0.5
        if self.nbd_site == 'c3':
            plt.plot(self.options.tspan, self.nbd_avgs[0], 'r.',
                     label='c3 data', alpha=alpha)
        elif self.nbd_site == 'c62':
            plt.plot(self.options.tspan, self.nbd_avgs[1], 'g.',
                     label='c62 data', alpha=alpha)
        #plt.plot(nbd.time_other, nbd_avgs[2], 'b.', label='c120 data',
        #          alpha=alpha)
        #plt.plot(nbd.time_other, nbd_avgs[3], 'm.', label='c122 data',
        #         alpha=alpha)
        #plt.plot(nbd.time_other, nbd_avgs[4], 'k.', label='c126 data',
        #         alpha=alpha)

    def get_observable_timecourse(self, position):
        """Gets the timecourse associated with the experimental observable."""
        # TODO Need to be able to get the indices for scaling parameters
        # from the model so that they're not hardcoded

        # Run the simulation and scale the simulated timecourse
        yout = self.simulate(position, observables=True)
        params = self.cur_params(position)
        timecourse = ((yout[self.nbd_observable] /
                       self.options.model.parameters['Bax_0'].value)
                      * params[3])
        return timecourse

    def fit_plotting_function(self, position):
        """Gets the observable timecourse and plots it against the data."""
        timecourse = self.get_observable_timecourse(position)

        # Make the plot
        plt.figure()
        self.plot_data()
        plt.plot(self.options.tspan, timecourse)
        plt.xlabel('Time')
        plt.ylabel('Concentration')
        fontP = FontProperties() 
        fontP.set_size('small')
        plt.legend(loc='upper center', prop=fontP, ncol=5,
                    bbox_to_anchor=(0.5, 1.1), fancybox=True, shadow=True)

    # A function to generate the likelihood function
    def get_likelihood_function(self):
        """Returns a likelihood function for the specified NBD site."""
        if self.nbd_site == 'c3':
            data_index = 0
        elif self.nbd_site == 'c62':
            data_index = 1
        else:
            raise Exception('Invalid value for nbd_site!')

        def likelihood(mcmc, position):
            yout = mcmc.simulate(position, observables=True)

            # TODO Need to be able to get the indices from the model so that
            # they're not hardcoded
            params = mcmc.cur_params(position)

            timecourse = ((yout[self.nbd_observable] /
                           mcmc.options.model.parameters['Bax_0'].value)
                          * params[3])

            return numpy.sum((self.nbd_avgs[data_index] - timecourse)**2 /
                             (2 * self.nbd_stds[data_index]**2))

        return likelihood

    def get_basename(self, model_name):
        """A function for standardizing, in one place, the format for pickled
        NBD_MCMC objects.
        """
        return '%s_%s_%s_%d_s%d' % (model_name, self.nbd_site,
                                    self.nbd_observable, self.options.nsteps,
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

        basename = '%s_%s_%s_%d_s%d.xxx' % (model, nbd_site, nbd_observable,
                                        nsteps, random_seed)

    With the suffix ``.xxx`` separated by a dot and the seed coming last in
    the underscore-separated arguments.

    Parameter
    ---------
    filenames : list of strings
        List of strings representing the chain filenames to be sorted into
        groups, e.g., of the type returned by ``glob.glob()``.

    Returns
    -------
    dict of lists of MCMC filenames. The keys in the dict are the filename
    prefixes that represent the arguments to the MCMC procedure (e.g.,
    ``tard_c3_iBax_4000``. Each dict entry contains a list of MCMC filenames
    associated with those run conditions.
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
